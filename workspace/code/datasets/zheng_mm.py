"""
Convert Zheng 2021 (GSE156728) pan-cancer T-cell counts → AnnData H5AD.

Input:
    shared/data/zheng2021/GSE156728_metadata.txt.gz
    shared/data/zheng2021/GSE156728_<CANCER>_10X.CD{4,8}.counts.txt.gz

Output:
    shared/data/zheng_mm_tcell.h5ad  (melanoma metastases, CD4+CD8 merged)

Rows=genes, cols=cells format. First row=header (cell IDs), first col=gene names.
Metadata has: cancerType, patient, libraryID, loc, meta.cluster, platform.

Rare populations of interest:
    CD8.c12.Tex.CXCL13  — TPEX / CXCR5+ precursor-exhausted CD8 T
    CD8.c11.Tex.PDCD1   — terminally exhausted CD8
    CD4.c20.Treg.TNFRSF9 — suppressive Treg
    CD8.c16.MAIT.SLC4A10 — MAIT
"""

from __future__ import annotations

import gzip
import pathlib

import numpy as np
import pandas as pd


def load_counts(path: pathlib.Path):
    """Read Zheng 2021 gene-by-cell gz text into (cell_ids, gene_names, X_cells_x_genes)."""
    with gzip.open(path, "rt") as fh:
        header = fh.readline().strip().split("\t")
    cell_ids = header  # cols are cell IDs
    print(f"  {path.name}: {len(cell_ids)} cells")

    # Load as DataFrame
    df = pd.read_csv(path, sep="\t", index_col=0)
    gene_names = list(df.index)
    X = df.to_numpy(dtype=np.float32).T  # cells × genes
    return cell_ids, gene_names, X


def main():
    src = pathlib.Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/data/zheng2021")

    print("Loading metadata…")
    meta = pd.read_csv(src / "GSE156728_metadata.txt.gz", sep="\t", index_col=0)
    print(f"  {meta.shape}")

    print("Loading MM CD8…")
    cid_cd8, genes_cd8, X_cd8 = load_counts(src / "GSE156728_MM_10X.CD8.counts.txt.gz")

    print("Loading MM CD4…")
    cid_cd4, genes_cd4, X_cd4 = load_counts(src / "GSE156728_MM_10X.CD4.counts.txt.gz")

    # Intersect genes
    inter = sorted(set(genes_cd8) & set(genes_cd4))
    print(f"  gene intersection: {len(inter)}")

    g8i = {g: i for i, g in enumerate(genes_cd8)}
    g4i = {g: i for i, g in enumerate(genes_cd4)}
    cols_cd8 = np.array([g8i[g] for g in inter], dtype=np.int64)
    cols_cd4 = np.array([g4i[g] for g in inter], dtype=np.int64)

    X_all = np.vstack([X_cd8[:, cols_cd8], X_cd4[:, cols_cd4]])
    all_cids = list(cid_cd8) + list(cid_cd4)
    print(f"  merged: {X_all.shape}")

    # Metadata reindex
    meta_aligned = meta.reindex(all_cids)
    missing = int(meta_aligned.isna().all(axis=1).sum())
    print(f"  metadata missing for {missing} cells; keeping them with NaN")

    import anndata as ad
    import scipy.sparse as sp
    a = ad.AnnData(X=sp.csr_matrix(X_all))
    a.obs_names = pd.Index(all_cids)
    a.var_names = pd.Index(inter)
    a.var["gene_symbol"] = a.var_names.astype(str).to_numpy()
    for col in meta_aligned.columns:
        a.obs[col] = meta_aligned[col].values
    a.obs["batch"] = a.obs["patient"].astype(str)
    a.obs["cell_type"] = a.obs["meta.cluster"].astype(str)
    a.obs["labels"] = a.obs["meta.cluster"].astype(str)

    out = pathlib.Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/data/zheng_mm_tcell.h5ad")
    a.write_h5ad(out)
    print(f"\nWrote {out}")
    print(f"  shape {a.shape}")
    print(f"  patient count: {a.obs['batch'].nunique()}")
    print(f"  meta.cluster top:")
    print(a.obs["meta.cluster"].value_counts().head(15).to_string())


if __name__ == "__main__":
    main()
