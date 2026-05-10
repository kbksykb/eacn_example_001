"""
Synthetic multi-batch scRNA-seq simulator for pipeline validation.

Generates a "pancreas-like" dataset: 4 cell types at 40/30/25/5% ratios
(alpha / beta / delta / epsilon), 3 batches with batch-specific shifts.
The 5% epsilon class is our label-holdout rare subpopulation.

Counts are NB-distributed per cell; genes are type-specific markers + shared
housekeeping + batch-specific perturbations. Resulting AnnData has:
    n_cells ≈ 6,000   (scale: 3 batches × ~2,000 cells)
    n_genes = 2,000
    obs['batch']           — str
    obs['cell_type_public']— str (will be masked for rare holdout)
    obs['heldout_rare']    — bool (True for epsilon)

This is ONLY for pipeline-plumbing validation. Real benchmark runs use the
scIB Baron/Muraro/Segerstolpe/Wang pancreas dataset once network lands.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class SimConfig:
    n_cells_per_batch: int = 2000
    n_batches: int = 3
    n_genes: int = 2000
    n_cell_types: int = 4
    abundances: tuple = (0.40, 0.30, 0.25, 0.05)  # 4th is epsilon
    rare_label: str = "epsilon"
    labels: tuple = ("alpha", "beta", "delta", "epsilon")
    marker_strength: float = 2.5
    batch_shift_sd: float = 0.5
    seed: int = 0


def simulate(cfg: SimConfig | None = None):
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp

    cfg = cfg or SimConfig()
    assert len(cfg.abundances) == cfg.n_cell_types
    assert len(cfg.labels) == cfg.n_cell_types
    rng = np.random.default_rng(cfg.seed)

    n_total = cfg.n_cells_per_batch * cfg.n_batches
    labels = []
    batches = []
    for b in range(cfg.n_batches):
        draw = rng.choice(cfg.n_cell_types, size=cfg.n_cells_per_batch, p=cfg.abundances)
        labels.append(draw)
        batches.append(np.full(cfg.n_cells_per_batch, b, dtype=np.int32))
    labels = np.concatenate(labels)
    batches = np.concatenate(batches)

    # Gene blocks: 100 markers per cell type + ambient housekeeping.
    # Epsilon (last type) gets a smaller marker block and moderate overlap with alpha
    # so it sits at the boundary of the abundant alpha population — the realistic
    # "absorbable rare" case.
    markers_per_type = 100
    marker_start = [i * markers_per_type for i in range(cfg.n_cell_types)]
    marker_end = [(i + 1) * markers_per_type for i in range(cfg.n_cell_types)]

    # Mean expression (log-space) — 0.5 baseline, +marker_strength for type markers
    base_mu = np.full((n_total, cfg.n_genes), 0.5, dtype=np.float32)
    for t in range(cfg.n_cell_types):
        mask = labels == t
        # Full marker strength on own block
        base_mu[np.ix_(mask, np.arange(marker_start[t], marker_end[t]))] += cfg.marker_strength
        # Epsilon (last type): also partially express alpha's markers at ~40% strength
        if t == cfg.n_cell_types - 1:
            base_mu[np.ix_(mask, np.arange(marker_start[0], marker_end[0]))] += 0.4 * cfg.marker_strength
            # And down-weight its own block slightly to make it harder
            base_mu[np.ix_(mask, np.arange(marker_start[t], marker_end[t]))] -= 0.5 * cfg.marker_strength

    # Batch-specific additive shift across shared housekeeping genes (not the markers)
    housekeeping = np.arange(cfg.n_cell_types * markers_per_type, cfg.n_genes)
    batch_shift = rng.normal(0, cfg.batch_shift_sd, size=(cfg.n_batches, housekeeping.size)).astype(np.float32)
    for b in range(cfg.n_batches):
        mask = batches == b
        base_mu[np.ix_(mask, housekeeping)] += batch_shift[b]

    # NB counts: X ~ Poisson(exp(mu))
    mu = np.exp(base_mu)
    X = rng.poisson(mu).astype(np.float32)

    obs = pd.DataFrame({
        "batch": pd.Categorical([f"batch_{b}" for b in batches]),
        "cell_type_public": pd.Categorical([cfg.labels[t] for t in labels]),
        "__hidden_label__": [cfg.labels[t] for t in labels],
        "heldout_rare": (labels == (cfg.n_cell_types - 1)),  # last cat = rare (epsilon)
        "oracle_rare_score": np.where(labels == (cfg.n_cell_types - 1), 2.0, 0.0).astype(np.float32),
    })
    obs_names = [f"cell_{i:06d}" for i in range(n_total)]
    var_names = [f"g_{j:04d}" for j in range(cfg.n_genes)]

    a = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var={"gene_ids": var_names})
    a.obs_names = obs_names
    a.var_names = var_names

    # Mask the rare label (keep cells, null the public label)
    new = a.obs["cell_type_public"].astype(object).copy()
    new[a.obs["heldout_rare"]] = np.nan
    a.obs["cell_type_public"] = new.astype("category")

    logger.info(
        "simulated %d cells × %d genes, %d batches, abundances=%s (rare=%s)",
        a.n_obs, a.n_vars, cfg.n_batches, cfg.abundances, cfg.rare_label,
    )
    return a


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    a = simulate()
    print(a)
    print("obs head:")
    print(a.obs.head())
    print("rare cells:", int(a.obs["heldout_rare"].sum()))
