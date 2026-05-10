"""
Convert Cheng 2021 PAAD .csv.gz → AnnData H5AD.

Source: GSE154763 (Cheng et al. 2021 Cell, pan-cancer myeloid atlas).
PAAD = pancreatic ductal adenocarcinoma slice.

Output: /mnt/.../shared/data/cheng_paad_myeloid.h5ad with
    .X = log-normalized expression
    .obs['patient_id'] (== batch for integration)
    .obs['cell_type'] (myeloid subtype from author annotation)
    .var['gene_symbol'] (gene symbols — essential for Immunology oracle scoring)
"""

from __future__ import annotations

import gzip
import pathlib

import numpy as np
import pandas as pd


def main():
    src = pathlib.Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/data/cheng2021")
    meta_path = src / "GSE154763_PAAD_metadata.csv.gz"
    expr_path = src / "GSE154763_PAAD_normalized_expression.csv.gz"
    out_path = src.parent / "cheng_paad_myeloid.h5ad"

    print(f"Loading metadata from {meta_path}…")
    meta = pd.read_csv(meta_path, index_col=0)
    print(f"  metadata shape {meta.shape}")
    print(f"  columns: {list(meta.columns)[:10]}")
    print(f"  head:\n{meta.head(3)}")

    print(f"Loading expression from {expr_path}…")
    expr = pd.read_csv(expr_path, index_col=0)
    print(f"  expr shape {expr.shape}")
    print(f"  first index: {expr.index[0]}, first column: {expr.columns[0]}")

    # Figure out orientation: Cheng stores as genes × cells
    # if rows are genes (13k) and cols are cells (7-8k)
    if expr.shape[0] > expr.shape[1]:
        # rows = genes
        gene_names = expr.index
        cell_ids = expr.columns
        X = expr.T.to_numpy(dtype=np.float32)
    else:
        gene_names = expr.columns
        cell_ids = expr.index
        X = expr.to_numpy(dtype=np.float32)
    print(f"  transposed to ({len(cell_ids)} cells × {len(gene_names)} genes)")

    # Align meta to cell_ids
    meta = meta.reindex([str(c) for c in cell_ids])
    print(f"  metadata aligned; meta.index NaN count: {meta.isna().all(axis=1).sum()}")

    import anndata as ad
    import scipy.sparse as sp
    a = ad.AnnData(X=sp.csr_matrix(X))
    a.obs_names = pd.Index([str(c) for c in cell_ids])
    a.var_names = pd.Index([str(g) for g in gene_names])
    a.var["gene_symbol"] = a.var_names.astype(str).to_numpy()
    for col in meta.columns:
        a.obs[col] = meta[col].values

    # Identify batch column — likely 'patient' or 'sample' or 'Patient'
    obs_cols = {c.lower(): c for c in a.obs.columns}
    for key in ("patient", "sample", "donor", "batch", "sampleid", "patient_id"):
        if key in obs_cols:
            a.obs["batch"] = a.obs[obs_cols[key]].astype(str)
            print(f"  batch key found: {obs_cols[key]}")
            break

    # Cell-type column
    for key in ("celltype", "cell_type", "subcluster", "cluster", "sub_celltype"):
        if key in obs_cols:
            a.obs["cell_type"] = a.obs[obs_cols[key]].astype(str)
            print(f"  cell_type key found: {obs_cols[key]}")
            break

    a.write_h5ad(out_path)
    print(f"Wrote {out_path}")
    print(f"  shape {a.shape}")
    print(f"  obs keys: {list(a.obs.columns)[:15]}")
    if "batch" in a.obs.columns:
        print(f"  batch counts (top): {a.obs['batch'].value_counts().head(5).to_dict()}")
    if "cell_type" in a.obs.columns:
        print(f"  cell_type counts: {a.obs['cell_type'].value_counts().head(10).to_dict()}")


if __name__ == "__main__":
    main()
