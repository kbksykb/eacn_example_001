"""scVI wrapper."""

from __future__ import annotations
import numpy as np


def run(adata, batch_key: str, seed: int) -> np.ndarray:
    import scvi
    import scanpy as sc

    a = adata.copy()
    sc.pp.highly_variable_genes(
        a, n_top_genes=2000, batch_key=batch_key, flavor="seurat_v3"
    )
    a = a[:, a.var["highly_variable"]].copy()

    scvi.settings.seed = seed
    scvi.model.SCVI.setup_anndata(a, batch_key=batch_key)
    model = scvi.model.SCVI(a, n_latent=30)
    model.train(max_epochs=200, early_stopping=True, check_val_every_n_epoch=10)
    return model.get_latent_representation()
