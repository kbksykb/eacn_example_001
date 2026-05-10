"""Mathematics-owned L_rare term for RareShield.

Mathematics agent (agent-mozusggu). 2026-05-10 (v2 bug-fix).
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

Mathematical form (per Math's colm_admissibility.md §3, two-sided |τ| variant):

    L_rare = E_w[ τ_stat(C_w)² ]

where τ(x) = log(m_post(N(T(x), r_N)) / m_pre(N(x, r_N))) and m_* is the Gaussian-kernel
local mass with a **fixed, shared** bandwidth r_N derived from z_pre. Using a fixed
(non-adaptive) bandwidth is essential: if bandwidth adapts to the local k-th-NN
distance in z_post, a tight rare cluster crushed into a denser region will still
register as "dense" because its k-NN distance remains small, yielding tau ≈ 0
despite the collapse. The fixed bandwidth anchors both measures to the same
physical scale so post mass reflects actual relocation (CompBio bug report,
2026-05-10).

The gated variant `τ² · ReLU(0.4 − L(w))` that used a saturating P_d proxy has been
replaced with `τ²` directly: |τ| ≈ 0 under admissible transport (no penalty regression)
and |τ| ≥ Ω(1) under collapse. The gate was redundant at the cost of a fragile proxy.

Usage:

    import torch
    from compute_l_rare import compute_l_rare, LRareConfig

    cfg = LRareConfig(k_N=50)
    loss_rare, diag = compute_l_rare(z_pre, z_post, motif_ids, cfg)
    total_loss = l_scvi + lambda_m * loss_rare + lambda_w * l_within + lambda_b * l_between
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch


@dataclass
class LRareConfig:
    k_N: int = 50                   # kNN size for global bandwidth determination
    bandwidth_quantile: float = 0.5  # use the median of pre k-NN distances as bandwidth
    eps: float = 1e-8               # numerical stability for log
    max_cells_per_motif: int = 2048 # sub-sample motifs larger than this


def _pairwise_sqdist(a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
    """Memory-light pairwise squared Euclidean distance."""
    return (a.unsqueeze(1) - b.unsqueeze(0)).pow(2).sum(-1)


def _global_bandwidth(z: torch.Tensor, k: int, quantile: float) -> torch.Tensor:
    """Quantile of k-th nearest-neighbor squared-distances across all cells.

    Returns a scalar h² used as a fixed Gaussian-kernel bandwidth so both pre- and
    post-integration local-mass counts are measured at the same physical scale.
    """
    with torch.no_grad():
        d2 = _pairwise_sqdist(z, z)
        kth, _ = torch.kthvalue(d2, k=k + 1, dim=1)  # +1 because self-distance at k=1
        return kth.quantile(quantile)


def _local_mass(z_anchor: torch.Tensor, z_cloud: torch.Tensor,
                h2: torch.Tensor) -> torch.Tensor:
    """Gaussian-kernel local mass around each row of z_anchor in z_cloud.

    Uses FIXED bandwidth h² rather than the adaptive k-th-NN-distance from before.
    Returns (n_anchor,) tensor; gradient flows through z_anchor if it requires grad.
    """
    d2 = _pairwise_sqdist(z_anchor, z_cloud)
    w = torch.exp(-d2 / (2.0 * h2))
    return w.sum(dim=1)


def _compute_tau(z_pre: torch.Tensor, z_post: torch.Tensor,
                 cfg: LRareConfig) -> torch.Tensor:
    """Per-cell τ(x) = log(m_post(N(T(x), r_N)) / m_pre(N(x, r_N))).

    Uses a FIXED bandwidth derived from z_pre so the two masses are measured at
    the same physical scale. z_pre is detached; gradient flows through z_post.
    """
    z_pre = z_pre.detach()
    h2 = _global_bandwidth(z_pre, cfg.k_N, cfg.bandwidth_quantile)

    # Pre-mass uses z_pre for both anchor and cloud (all detached) — scalar constant
    m_pre = _local_mass(z_pre, z_pre, h2).detach()

    # Post-mass: anchor is z_post (grad-enabled), cloud is z_post (detached)
    # The detach on the cloud makes the kernel position-shift cleanly attributable
    # to the anchor coordinates; this is the standard practice for density-based
    # losses (cf. Neural OT, Sinkhorn with stop-grad on targets).
    m_post = _local_mass(z_post, z_post.detach(), h2)

    return torch.log(m_post + cfg.eps) - torch.log(m_pre + cfg.eps)


def compute_l_rare(
    z_pre: torch.Tensor,
    z_post: torch.Tensor,
    motif_ids: torch.Tensor,
    cfg: Optional[LRareConfig] = None,
) -> tuple[torch.Tensor, dict]:
    """L_rare = E_w[ τ_stat(C_w)² ] over motifs.

    Arguments
    ---------
    z_pre : (n, d) — pre-integration embedding, detached internally.
    z_post : (n, d) — post-integration embedding; grad flows through.
    motif_ids : (n,) long — motif assignment per cell, -1 for non-motif.
    cfg : LRareConfig.

    Returns
    -------
    loss : scalar tensor (0 if no motifs).
    diagnostics : per-motif dict with tau_mean, tau_abs_mean, term, n_cells.
    """
    cfg = cfg or LRareConfig()
    assert z_pre.shape == z_post.shape, "embeddings must have the same shape"
    assert motif_ids.shape[0] == z_pre.shape[0], "motif_ids must align with embeddings"

    tau_all = _compute_tau(z_pre, z_post, cfg)

    device = z_post.device
    unique_motifs = motif_ids.unique()
    unique_motifs = unique_motifs[unique_motifs >= 0]  # drop -1 (non-motif)

    if unique_motifs.numel() == 0:
        return torch.zeros((), device=device, requires_grad=True), {"n_motifs": 0}

    loss_terms = []
    diag = {"per_motif": []}
    for w_id in unique_motifs.tolist():
        mask = (motif_ids == w_id)
        if mask.sum() == 0:
            continue
        tau_w = tau_all[mask]
        if tau_w.numel() > cfg.max_cells_per_motif:
            idx = torch.randperm(tau_w.numel(), device=device)[: cfg.max_cells_per_motif]
            tau_w = tau_w[idx]

        # Two-sided |τ| handles both absorption (τ>0) and dispersion (τ<0).
        # Trimmed mean is robust to outliers; squaring gives the quadratic penalty.
        tau_sorted, _ = torch.sort(tau_w)
        nq = tau_sorted.numel()
        q_lo = tau_sorted[int(0.1 * nq):int(0.9 * nq) + 1]
        tau_stat_sq = q_lo.pow(2).mean()

        loss_terms.append(tau_stat_sq)
        diag["per_motif"].append({
            "motif_id": w_id,
            "n_cells": int(mask.sum()),
            "tau_mean": float(tau_w.mean().item()),
            "tau_abs_mean": float(tau_w.abs().mean().item()),
            "term": float(tau_stat_sq.item()),
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
    motif_ids = torch.full((n,), -1, dtype=torch.long)
    rare = torch.randperm(n)[:40]
    motif_ids[rare] = 0
    # Rare cells offset by 3σ in dim 0 from the abundant cloud
    z_pre[rare] = z_pre[rare] + torch.tensor([3.0] + [0.0] * (d - 1))

    cfg = LRareConfig(k_N=50)

    # Case A: well-preserved
    z_post_good = (z_pre + 0.1 * torch.randn(n, d)).detach().requires_grad_()
    loss_good, diag_good = compute_l_rare(z_pre, z_post_good, motif_ids, cfg)
    print(f"L_rare (preserved):  {loss_good.item():.4f}   "
          f"(|tau| mean = {diag_good['per_motif'][0]['tau_abs_mean']:.3f})")

    # Case B: crushed — rare cells moved onto abundant centroid
    z_post_bad = z_pre.clone()
    z_post_bad[rare] = 0.3 * torch.randn(len(rare), d)
    z_post_bad = z_post_bad.detach().requires_grad_()
    loss_bad, diag_bad = compute_l_rare(z_pre, z_post_bad, motif_ids, cfg)
    print(f"L_rare (crushed):    {loss_bad.item():.4f}   "
          f"(|tau| mean = {diag_bad['per_motif'][0]['tau_abs_mean']:.3f})")

    loss_bad.backward()
    print(f"grad on crushed motif (mean |grad|): "
          f"{z_post_bad.grad[rare].abs().mean().item():.4e}")

    assert loss_bad.item() > 10 * loss_good.item(), \
        "Self-test failed: crushed should be >> preserved"
    print("PASS: crushed >> preserved, gradient flows through rare cells.")
