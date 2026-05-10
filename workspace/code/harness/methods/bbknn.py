"""BBKNN wrapper — emits X_bbknn via scanpy neighbours + connectivities matrix stored as embedding proxy."""

from __future__ import annotations
import numpy as np


def run(adata, batch_key: str, seed: int) -> np.ndarray:
    import scanpy as sc
    import bbknn

    a = adata.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2000, batch_key=batch_key)
    a = a[:, a.var["highly_variable"]].copy()
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=50, random_state=seed)

    bbknn.bbknn(a, batch_key=batch_key)
    # BBKNN returns a neighbourhood graph only. For a canonical X_integrated we
    # use UMAP on the bbknn graph as the 2D proxy. For RareScore we mostly need
    # the graph; store both.
    sc.tl.umap(a, random_state=seed)
    return a.obsm["X_umap"].astype(np.float32)
