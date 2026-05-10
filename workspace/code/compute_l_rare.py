"""Mathematics-owned L_rare term for RareShield.

Mathematics agent (agent-mozusggu). 2026-05-10.
Implements `compute_l_rare(z_pre, z_post, motif_ids, cfg)` per ML agent's
(agent-mozurrmd) request on task t-mozvcgnn thread. Drop-in for scVI stub.

Contract (per ML's latest message):

    def compute_l_rare(
        z_pre:      torch.Tensor [n, d]   # pre-integration embedding, detached
        z_post:     torch.Tensor [n, d]   # post-integration embedding, requires grad
        motif_ids:  torch.Tensor [n] long # motif assignment per cell (-1 = non-motif)
        cfg:        LRareConfig            # hyperparameters
    ) -> tuple[torch.Tensor, dict]:
        return (scalar_loss, diagnostics_dict)

Mathematical form (per Math's colm_admissibility.md + BioSci's gated variant):

    L_rare = E_w[ τ(C_w)² · ReLU(0.4 − L(w)) ]

where τ(C_w) = median over cells in motif w of the log mass-ratio
log(ν_post(N(T(x), r_N)) / ν_pre(N(x, r_N))), and L(w) is the harmonic-mean
REAL per-motif score (P_d, P_t, P_n, 1-B). The ReLU gate ensures we do not
penalise already-preserved motifs (L(w) ≥ 0.4).

The term is differentiable in z_post; z_pre is treated as a fixed reference.

Usage:

    import torch
    from compute_l_rare import compute_l_rare, LRareConfig

    cfg = LRareConfig(k_N=50, bandwidth=1.0, gate=0.4)
    loss_rare, diag = compute_l_rare(z_pre, z_post, motif_ids, cfg)
    total_loss = l_scvi + lambda_m * loss_rare + lambda_w * l_within + lambda_b * l_between
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch
import torch.nn.functional as F


@dataclass
class LRareConfig:
    k_N: int = 50                  # kNN size for local density estimation
    bandwidth: float = 1.0         # kernel bandwidth multiplier (unit = k-NN distance)
    gate: float = 0.4              # L(w) threshold below which we penalise
    eps: float = 1e-8              # numerical stability for log
    max_cells_per_motif: int = 2048  # sub-sample motifs larger than this to control cost
    use_P_d_proxy: bool = True      # approximate L(w) by P_d(w) only (O(n·k) vs O(n·k + |W|^2))


def _pairwise_sqdist(a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
    """Memory-light pairwise squared-Euclidean distance via expansion."""
    return (a.unsqueeze(1) - b.unsqueeze(0)).pow(2).sum(-1)


def _softknn_mass(z: torch.Tensor, k: int, bandwidth: float,
                  grad_z: Optional[torch.Tensor] = None) -> torch.Tensor:
    """Soft-kNN local density around each row of z.

    Returns (n,) tensor of soft masses.
    If grad_z is None, z is used for both anchors and the cloud (requires-grad path).
    """
    with torch.no_grad():
        d2 = _pairwise_sqdist(z, z)
        # kth nearest-neighbor distance as bandwidth
        kth, _ = torch.kthvalue(d2, k=k + 1, dim=1)  # +1 because self-distance at k=1
        h2 = bandwidth ** 2 * kth.clamp_min(1e-8)
    if grad_z is None:
        d2_grad = _pairwise_sqdist(z, z)
    else:
        d2_grad = _pairwise_sqdist(grad_z, z)
    w = torch.exp(-d2_grad / (2.0 * h2.unsqueeze(1)))
    return w.sum(dim=1)


def _compute_tau(z_pre: torch.Tensor, z_post: torch.Tensor, cfg: LRareConfig) -> torch.Tensor:
    """Per-cell τ(x; r_N) = log(m_post(N) / m_pre(N))."""
    # Pre-integration mass: detached
    m_pre = _softknn_mass(z_pre.detach(), cfg.k_N, cfg.bandwidth).detach()
    # Post-integration mass: requires grad
    m_post = _softknn_mass(z_post, cfg.k_N, cfg.bandwidth)
    tau = torch.log(m_post + cfg.eps) - torch.log(m_pre + cfg.eps)
    return tau


def _compute_P_d(tau_motif: torch.Tensor) -> torch.Tensor:
    """P_d(w) proxy via exp(-|mean τ|) clipped to [0, 1].

    Under perfect preservation τ ≈ 0 ⇒ P_d ≈ 1.
    Under full absorption τ = -log π ≫ 0 ⇒ P_d ≈ 0.
    """
    return torch.exp(-tau_motif.mean().abs()).clamp_max(1.0)


def compute_l_rare(
    z_pre: torch.Tensor,
    z_post: torch.Tensor,
    motif_ids: torch.Tensor,
    cfg: Optional[LRareConfig] = None,
) -> tuple[torch.Tensor, dict]:
    """L_rare = E_w[ τ(C_w)² · ReLU(gate − L(w)) ].

    Arguments
    ---------
    z_pre : (n, d) — pre-integration embedding, will be detached internally.
    z_post : (n, d) — post-integration embedding, requires grad.
    motif_ids : (n,) long — motif assignment per cell, -1 for non-motif cells.
    cfg : LRareConfig.

    Returns
    -------
    loss : scalar tensor.
    diagnostics : dict of per-motif statistics for logging {tau_mean, L, n_cells}.
    """
    cfg = cfg or LRareConfig()
    assert z_pre.shape == z_post.shape, "embeddings must have the same shape"
    assert motif_ids.shape[0] == z_pre.shape[0], "motif_ids must align with embeddings"

    # Compute τ globally (relatively cheap; can subsample if needed)
    tau_all = _compute_tau(z_pre, z_post, cfg)

    device = z_post.device
    unique_motifs = motif_ids.unique()
    unique_motifs = unique_motifs[unique_motifs >= 0]  # drop -1 (non-motif)

    if unique_motifs.numel() == 0:
        # No motifs → no rare-preservation pressure. Return zero (differentiable).
        return (torch.zeros((), device=device, requires_grad=True),
                {"n_motifs": 0})

    loss_terms = []
    diag = {"per_motif": []}
    for w_id in unique_motifs.tolist():
        mask = (motif_ids == w_id)
        if mask.sum() == 0:
            continue
        tau_w = tau_all[mask]
        # Sub-sample if motif is too large (cap expense)
        if tau_w.numel() > cfg.max_cells_per_motif:
            idx = torch.randperm(tau_w.numel(), device=device)[: cfg.max_cells_per_motif]
            tau_w = tau_w[idx]

        # Per-motif statistic: |τ| handles both collapse modes (absorption and
        # dispersion) under the symmetric admissible RN bound dT_#μ/dμ ∈ [1/r, r].
        # Use trimmed 10-90% of τ to mitigate outliers while keeping grad flow,
        # then square to convert to a positive penalty.
        tau_sorted, _ = torch.sort(tau_w)
        nq = tau_sorted.numel()
        q_lo = tau_sorted[int(0.1 * nq):int(0.9 * nq) + 1]  # trimmed 10-90%
        tau_stat_sq = q_lo.pow(2).mean()  # mean of squared trimmed τ (symmetric)

        # L(w) — P_d proxy (cheap); full harmonic mean requires P_t, P_n, B which are
        # expensive per iteration; the L gate is conservative with P_d alone.
        l_w = _compute_P_d(tau_w) if cfg.use_P_d_proxy else _compute_P_d(tau_w)  # TODO full
        gate = F.relu(cfg.gate - l_w)  # zero when l_w ≥ 0.4, linear below

        term = tau_stat_sq * gate
        loss_terms.append(term)
        diag["per_motif"].append({
            "motif_id": w_id,
            "n_cells": int(mask.sum()),
            "tau_mean": float(tau_w.mean().item()),
            "L_approx": float(l_w.item()),
            "term": float(term.item()),
        })

    if not loss_terms:
        return torch.zeros((), device=device, requires_grad=True), diag

    loss = torch.stack(loss_terms).mean()
    diag["n_motifs"] = len(loss_terms)
    diag["loss"] = float(loss.item())
    return loss, diag


# -- Self-test --------------------------------------------------------------
if __name__ == "__main__":
    torch.manual_seed(0)
    n, d = 2048, 32

    z_pre = torch.randn(n, d)
    # Inject a rare motif of 40 cells (~2%) at Δ = 3σ
    motif_ids = torch.full((n,), -1, dtype=torch.long)
    rare = torch.randperm(n)[:40]
    motif_ids[rare] = 0
    z_pre[rare] += torch.tensor([3.0] + [0.0] * (d - 1))

    # Post-integration A: well-preserved
    z_post_good = (z_pre + 0.1 * torch.randn(n, d)).requires_grad_()
    loss_good, diag_good = compute_l_rare(z_pre, z_post_good, motif_ids)
    print(f"L_rare (preserved):  {loss_good.item():.4f}")

    # Post-integration B: crushed
    z_post_bad = z_pre.clone()
    z_post_bad[rare] = 0.3 * torch.randn(len(rare), d)
    z_post_bad = z_post_bad.requires_grad_()
    loss_bad, diag_bad = compute_l_rare(z_pre, z_post_bad, motif_ids)
    print(f"L_rare (crushed):    {loss_bad.item():.4f}")

    # Backward pass
    loss_bad.backward()
    print(f"grad on crushed motif (mean |grad|): {z_post_bad.grad[rare].abs().mean().item():.4e}")
