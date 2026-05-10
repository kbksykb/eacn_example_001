"""
Merge Cheng 2021 GSE154763 cancer-type slices into ONE pan-cancer myeloid atlas.

Input: 6 cancer-type pairs (metadata + normalized_expression) at shared/data/cheng2021/.
Output: /mnt/.../shared/data/cheng_pancancer_myeloid.h5ad

Batch key: 'batch' = patient (consistent with Cheng's original).
Cancer-type key: 'cancer' = {ESCA, KIDNEY, LYM, OV-FTC, PAAD, THCA}.
"""

from __future__ import annotations

import gzip
import pathlib

import numpy as np
import pandas as pd


def load_one(cancer: str):
    src = pathlib.Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/data/cheng2021")
    meta = pd.read_csv(src / f"GSE154763_{cancer}_metadata.csv.gz", index_col=0)
    expr = pd.read_csv(src / f"GSE154763_{cancer}_normalized_expression.csv.gz", index_col=0)
    print(f"  [{cancer}] meta={meta.shape} expr={expr.shape}")
    # Cheng format: rows=cells, cols=genes (always). X is already cells × genes.
    gene_names = list(expr.columns)
    cell_ids = list(expr.index)
    X = expr.to_numpy(dtype=np.float32)
    meta = meta.reindex([str(c) for c in cell_ids])
    return X, gene_names, cell_ids, meta


def main():
    cancers = ["PAAD", "THCA", "ESCA", "KIDNEY", "LYM", "OV-FTC"]
    Xs = []
    metas = []
    cell_idss = []
    all_genes = None

    for c in cancers:
        X, gnames, cids, m = load_one(c)
        m["cancer_type"] = c
        Xs.append((X, gnames))
        metas.append(m)
        cell_idss.append([f"{c}_{cid}" for cid in cids])

    # Gene union, then restrict to intersection for cleaner concat
    gene_sets = [set(gn) for X, gn in Xs]
    inter = sorted(set.intersection(*gene_sets))
    print(f"Gene intersection across {len(cancers)} cancers: {len(inter)}")

    # Build index-to-position map per cancer, aligned to inter
    import scipy.sparse as sp
    pieces = []
    for (X, gn), m, cids in zip(Xs, metas, cell_idss):
        gn_to_i = {g: i for i, g in enumerate(gn)}
        cols = np.array([gn_to_i[g] for g in inter], dtype=np.int64)
        X_sub = X[:, cols]
        pieces.append((X_sub, m))

    X_all = np.vstack([p[0] for p in pieces])
    meta_all = pd.concat([p[1] for p in pieces], axis=0, ignore_index=False)
    # Flatten cell IDs with cancer prefix to avoid collisions
    obs_names = []
    for cancer, cids in zip(cancers, cell_idss):
        obs_names.extend(cids)

    import anndata as ad
    a = ad.AnnData(X=sp.csr_matrix(X_all))
    a.obs_names = pd.Index(obs_names)
    a.var_names = pd.Index(inter)
    a.var["gene_symbol"] = a.var_names.astype(str).to_numpy()
    for col in meta_all.columns:
        a.obs[col] = meta_all[col].values

    # Canonicalize columns
    if "patient" in a.obs.columns:
        a.obs["batch"] = a.obs["patient"].astype(str)
    if "MajorCluster" in a.obs.columns:
        a.obs["cell_type"] = a.obs["MajorCluster"].astype(str)
        a.obs["labels"] = a.obs["MajorCluster"].astype(str)

    out = pathlib.Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/data/cheng_pancancer_myeloid.h5ad")
    a.write_h5ad(out)
    print(f"\nWrote {out}")
    print(f"  shape={a.shape}")
    print(f"  cancer_type: {a.obs['cancer_type'].value_counts().to_dict()}")
    print(f"  batch (patient) count: {a.obs['batch'].nunique()}")
    print(f"  MajorCluster: {a.obs['MajorCluster'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
