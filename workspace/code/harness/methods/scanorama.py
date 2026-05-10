"""Scanorama wrapper."""

from __future__ import annotations
import numpy as np


def run(adata, batch_key: str, seed: int) -> np.ndarray:
    import scanpy as sc
    import scanorama

    a = adata.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2000, batch_key=batch_key)
    a = a[:, a.var["highly_variable"]].copy()

    batches = a.obs[batch_key].astype(str)
    uniq = sorted(batches.unique())
    splits = [a[batches == b].copy() for b in uniq]
    corrected, _ = scanorama.correct_scanpy(splits, return_dimred=True)

    order = np.concatenate([np.where(batches == b)[0] for b in uniq])
    emb = np.empty((a.n_obs, corrected[0].obsm["X_scanorama"].shape[1]), dtype=np.float32)
    for b, piece in zip(uniq, corrected):
        emb[np.where(batches == b)[0]] = piece.obsm["X_scanorama"]
    return emb
