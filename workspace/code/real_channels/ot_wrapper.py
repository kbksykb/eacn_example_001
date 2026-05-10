"""
OT channel wrapper — calls ML's score_ot() with a two-sided test and surface the
p-value and statistic in the column names used by DS's Fig 2 rendering.

Until ML publishes a patched smoke test, we use the absolute τ for p-values, which
is strictly stronger: rare-cluster loss can manifest as either "mass destruction
into the void" (τ<<0 in ML's sign convention) or "absorption into dense regions"
(τ>>0). Detector should catch both. See message thread with agent-mozurrmd.
"""

from __future__ import annotations

import numpy as np

try:
    from . import ot_channel as _ot  # type: ignore
    OTChannelConfig = _ot.OTChannelConfig
    _ml_score = _ot.score_ot
    OT_AVAILABLE = True
except ImportError:  # pragma: no cover
    OTChannelConfig = None
    _ml_score = None
    OT_AVAILABLE = False


def score_two_sided(emb_pre, emb_post, candidate_labels, batch_labels, cfg=None):
    """Two-sided variant of ML's score_ot: detects both mass-destruction and
    absorption. Returns the same dict shape but with p-values on |τ|."""
    if not OT_AVAILABLE:
        return {
            "channel": "ot_density_flow",
            "candidate_ids": np.array([], dtype=int),
            "stat": np.array([]),
            "p_values": np.array([]),
            "ci_low": np.array([]),
            "ci_high": np.array([]),
            "config": {},
        }

    # Call ML's module to get τ_obs + null
    # Easiest: patch the module-level routine by reusing its primitives
    import torch
    import numpy as _np

    cfg = cfg or OTChannelConfig()
    torch.manual_seed(cfg.seed)
    rng = _np.random.default_rng(cfg.seed)

    emb_pre_t = torch.as_tensor(emb_pre, dtype=cfg.dtype, device=cfg.device)
    emb_post_t = torch.as_tensor(emb_post, dtype=cfg.dtype, device=cfg.device)

    cand_ids = _np.unique(candidate_labels[candidate_labels >= 0])

    tau = _ot.per_cell_tau(emb_pre_t, emb_post_t, cfg.k_neighbors).cpu().numpy()
    abs_tau = _np.abs(tau)

    stat_obs = _np.empty(len(cand_ids))
    for i, cid in enumerate(cand_ids):
        mask = candidate_labels == cid
        stat_obs[i] = float(_np.median(abs_tau[mask]))

    stat_null = _np.empty((cfg.n_permutations, len(cand_ids)))
    for r in range(cfg.n_permutations):
        abs_tau_perm = rng.permutation(abs_tau)
        for i, cid in enumerate(cand_ids):
            mask = candidate_labels == cid
            stat_null[r, i] = float(_np.median(abs_tau_perm[mask]))

    p_values = (1 + _np.sum(stat_null >= stat_obs[None, :], axis=0)) / (cfg.n_permutations + 1)

    # Bootstrap CI on observed stat
    ci_low = _np.empty(len(cand_ids))
    ci_high = _np.empty(len(cand_ids))
    for i, cid in enumerate(cand_ids):
        mask = candidate_labels == cid
        tau_cand = abs_tau[mask]
        n_c = tau_cand.size
        if n_c < 2:
            ci_low[i], ci_high[i] = _np.nan, _np.nan
            continue
        boot = rng.choice(tau_cand, size=(200, n_c), replace=True)
        meds = _np.median(boot, axis=-1)
        ci_low[i], ci_high[i] = _np.quantile(meds, [0.025, 0.975])

    return {
        "channel": "ot_density_flow_two_sided",
        "candidate_ids": cand_ids,
        "stat": stat_obs,
        "signed_stat_mean": _np.array([
            float(_np.mean(tau[candidate_labels == cid])) for cid in cand_ids
        ]),  # the signed mean τ helps interpret (positive → absorbed, negative → scattered)
        "p_values": p_values,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "config": cfg.__dict__,
    }
