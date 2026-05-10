"""scANVI wrapper (semi-supervised)."""

from __future__ import annotations
import numpy as np


def run(adata, batch_key: str, seed: int, label_key: str = "cell_type") -> np.ndarray:
    import scvi
    import scanpy as sc

    a = adata.copy()
    sc.pp.highly_variable_genes(a, n_top_genes=2000, batch_key=batch_key, flavor="seurat_v3")
    a = a[:, a.var["highly_variable"]].copy()

    scvi.settings.seed = seed
    scvi.model.SCVI.setup_anndata(a, batch_key=batch_key, labels_key=label_key)
    svi = scvi.model.SCVI(a, n_latent=30)
    svi.train(max_epochs=100, check_val_every_n_epoch=10, early_stopping=True)

    scanvi = scvi.model.SCANVI.from_scvi_model(svi, unlabeled_category="Unknown", labels_key=label_key)
    scanvi.train(max_epochs=100, n_samples_per_label=100, check_val_every_n_epoch=10, early_stopping=True)
    return scanvi.get_latent_representation()
