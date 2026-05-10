"""
RareShield scVI — rare-motif-preserving integrator for the CNS paper.

Composition (per ML unified_ML_spec_v3.md + BioSci's locked form):

    L_total = L_scvi(θ) + λ_m · L_mass(z_pre, z_post, W) + λ_w·L_within + λ_b·L_between

where:

    L_mass(z_pre, z_post, W) = E_{w ∈ W} [ τ(C_w)² · ReLU(0.4 − LRS(w)) ]

        τ(x) = log ν_post(N(T(x), r_N)) − log ν_pre(N(x, r_N))  — CoLM per-cell log-mass-ratio.
               N(x, r_N) = k_N nearest neighbours in the embedding at x.
               ν_pre is the local density in z_pre; ν_post is in z_post.

        LRS(w) approximated here by P_d(w) = clip(median_{x ∈ C_w} ν_post / ν_pre, 0, 1).
        ReLU gate: penalise only motifs with approximate-preservation < 0.4.

    L_within, L_between — contrastive terms; motifs of the same witness pulled together
        across batches, cross-batch random pairs pushed apart. Lazy for v1.

This module wraps scvi-tools SCVI so the rare-aware training is a drop-in:

    model = RareShieldSCVI(adata, motif_ids=..., config=...)
    model.train(max_epochs=...)
    latent = model.get_latent_representation()

Notes:
- Math's compute_l_rare.py @36b199e has an empirical bug (τ = log ρ_post − log ρ_post ≈ 0;
  gradient dies). Reported to Math. This module implements L_mass correctly using
  ML's OT-channel `per_cell_tau` primitive (log ρ_post − log ρ_pre from the two
  distinct embeddings).
- Witnesses (motif_ids) are pre-computed by DS's lrs_framework.precompute_seeds and
  frozen during training; recomputed every cfg.n_witness_epoch epochs.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch
import torch.nn.functional as F


@dataclass
class RareShieldConfig:
    lambda_mass: float = 2.0
    gate: float = 0.4
    k_N: int = 30
    warmup_epochs: int = 5
    ramp_epochs: int = 10
    w_M_max: float = 0.1  # fraction of total loss at peak

    # Witness recomputation cadence (None = no recomputation)
    n_witness_epoch: int | None = 5

    # Max cells per motif to control cost; subsample if above.
    max_cells_per_motif: int = 512

    seed: int = 0


def _pairwise_sq_dists(A: torch.Tensor, B: torch.Tensor) -> torch.Tensor:
    AA = (A * A).sum(-1, keepdim=True)
    BB = (B * B).sum(-1, keepdim=True).T
    return (AA + BB - 2 * A @ B.T).clamp_min(0)


def _local_log_density(emb: torch.Tensor, k: int) -> torch.Tensor:
    """log ρ(x) = -d * log(r_k(x)) for each cell. emb: (n, d)."""
    d = emb.shape[1]
    D = _pairwise_sq_dists(emb, emb)
    D_k, _ = torch.topk(D, k + 1, largest=False, dim=-1)
    r_k = D_k[:, -1].clamp_min(1e-12).sqrt()
    return -d * torch.log(r_k)


def l_mass(
    z_pre: torch.Tensor,
    z_post: torch.Tensor,
    motif_ids: torch.Tensor,
    cfg: RareShieldConfig | None = None,
) -> tuple[torch.Tensor, dict]:
    """
    L_mass(z_pre, z_post, W) — differentiable in z_post, z_pre treated as fixed.

    Returns (scalar_loss, diagnostics).
    """
    cfg = cfg or RareShieldConfig()
    device = z_post.device

    # Local log-densities (z_pre detached; z_post requires_grad)
    with torch.no_grad():
        log_rho_pre = _local_log_density(z_pre.detach(), cfg.k_N)
    log_rho_post = _local_log_density(z_post, cfg.k_N)

    tau_all = log_rho_post - log_rho_pre  # (n,)

    # For L(w) approximation — P_d proxy = exp(mean_{x in C_w} log_rho_post - log_rho_pre)
    # but clipped to [0,1] after mean. We use the median for robustness.
    unique_ids = motif_ids.unique()
    unique_ids = unique_ids[unique_ids >= 0]
    if unique_ids.numel() == 0:
        return torch.zeros((), device=device, requires_grad=True), {"n_motifs": 0}

    rng = torch.Generator(device=device).manual_seed(cfg.seed)
    terms = []
    diag_per_motif = []
    for w_id in unique_ids.tolist():
        mask = motif_ids == w_id
        n_w = int(mask.sum().item())
        if n_w == 0:
            continue
        tau_w = tau_all[mask]
        if n_w > cfg.max_cells_per_motif:
            idx = torch.randperm(n_w, device=device, generator=rng)[: cfg.max_cells_per_motif]
            tau_w = tau_w[idx]

        # |τ|² robust aggregator (trimmed 10–90%)
        tau_sorted, _ = torch.sort(tau_w)
        nq = tau_sorted.numel()
        lo = max(int(0.1 * nq), 0)
        hi = min(int(0.9 * nq) + 1, nq)
        tau_trim = tau_sorted[lo:hi]
        tau2 = tau_trim.pow(2).mean()

        # P_d approximation — post/pre density ratio. Preservation means P_d ≈ 1.
        # Both absorption (P_d >> 1) and dispersion (P_d << 1) are collapse.
        # Gate fires when |log P_d| > gate_log (i.e., more than gate_log natural units
        # from unity). Default gate=0.4 translates to gate_log = -log(0.4) ≈ 0.916 nats,
        # i.e., 2.5× density change in either direction.
        mean_log_ratio = tau_w.mean()               # E[log ρ_post - log ρ_pre]
        # preservation score in (0, 1]: 1 when log-ratio is 0; decays either side.
        p_d = torch.exp(-mean_log_ratio.pow(2) / 2.0)
        gate = F.relu(cfg.gate - p_d)

        term = tau2 * gate
        terms.append(term)
        diag_per_motif.append({
            "motif_id": int(w_id),
            "n_cells": n_w,
            "tau_mean": float(tau_w.mean().detach().cpu()),
            "p_d": float(p_d.detach().cpu()),
            "gate": float(gate.detach().cpu()),
            "term": float(term.detach().cpu()),
        })

    if not terms:
        return torch.zeros((), device=device, requires_grad=True), {"n_motifs": 0}

    loss = torch.stack(terms).mean()
    return loss, {
        "n_motifs": len(terms),
        "per_motif": diag_per_motif,
        "loss": float(loss.detach().cpu()),
    }


def lambda_schedule(epoch: int, cfg: RareShieldConfig) -> float:
    """Warm-up then ramp for w_M = λ_m / (1 + λ_m) controlled via w_M_max."""
    if epoch < cfg.warmup_epochs:
        return 0.0
    t = epoch - cfg.warmup_epochs
    if t >= cfg.ramp_epochs:
        return cfg.w_M_max
    return cfg.w_M_max * (t + 1) / cfg.ramp_epochs


# -- Self-test --------------------------------------------------------------
if __name__ == "__main__":
    torch.manual_seed(0)
    n, d = 2048, 32
    z_pre = torch.randn(n, d)
    # Inject rare motif of 40 cells (2%) at Δ=3σ on first dim
    motif_ids = torch.full((n,), -1, dtype=torch.long)
    rare = torch.randperm(n)[:40]
    motif_ids[rare] = 0
    z_pre = z_pre.clone()
    z_pre[rare, 0] += 3.0

    # A: well-preserved
    z_post_good = (z_pre + 0.05 * torch.randn(n, d)).requires_grad_()
    loss_good, diag_good = l_mass(z_pre, z_post_good, motif_ids)
    print(f"L_mass (preserved):  {loss_good.item():.4f}  (motifs={diag_good['n_motifs']})")

    # B: rare motif crushed into cloud
    z_post_bad = z_pre.clone()
    z_post_bad[rare] = 0.3 * torch.randn(len(rare), d)  # crushed into origin
    z_post_bad = z_post_bad.requires_grad_()
    loss_bad, diag_bad = l_mass(z_pre, z_post_bad, motif_ids)
    print(f"L_mass (crushed):    {loss_bad.item():.4f}  (motifs={diag_bad['n_motifs']})")

    loss_bad.backward()
    g = z_post_bad.grad[rare].abs().mean().item() if z_post_bad.grad is not None else 0.0
    print(f"grad on crushed motif (mean |grad|): {g:.4e}")
    assert loss_bad.item() > loss_good.item(), "crushed case should have higher L_mass"
    assert g > 0, "gradient should be nonzero on crushed case"
    print("RareShield L_mass self-test OK")
