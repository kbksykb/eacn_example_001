"""Harmony (harmonypy) wrapper."""

from __future__ import annotations
import numpy as np


def run(adata, batch_key: str, seed: int) -> np.ndarray:
    import scanpy as sc
    import harmonypy as hm

    rng = np.random.default_rng(seed)
    _ = rng.bytes(1)  # keep rng reproducible usage

    a = adata.copy()
    if "X_pca" not in a.obsm:
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        sc.pp.highly_variable_genes(a, n_top_genes=2000, batch_key=batch_key)
        a = a[:, a.var["highly_variable"]].copy()
        sc.pp.scale(a, max_value=10)
        sc.tl.pca(a, n_comps=50, random_state=seed)

    ho = hm.run_harmony(a.obsm["X_pca"], a.obs, vars_use=[batch_key], random_state=seed)
    return np.asarray(ho.Z_corr.T)
