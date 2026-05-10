"""
Emit scIB-style metrics JSON per (method, dataset) for DS's Fig 3b.

Per DS's DATA_SPEC: shared/scib/<method>_<dataset>.json with keys
    {ari, nmi, asw, graph_conn, ilisi, clisi}

Quick implementation without full scib-metrics package (which has scikit-learn version
incompatibilities) — we compute each metric by hand using sklearn + scanpy.

For each integrated H5AD:
    ari, nmi : between Leiden clustering on X_integrated and true labels.
    asw      : silhouette score on X_integrated using true labels.
    graph_conn : graph connectivity (same labels, same-label-cells must stay connected).
    ilisi    : inverse LISI (integration LISI — higher is better at batch mixing).
    clisi    : cell-type LISI — lower is better at label separation.
"""

from __future__ import annotations

import json
import logging
import pathlib

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

SHARED = pathlib.Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared")


def compute_scib(adata, label_key: str = "labels", batch_key: str = "batch", seed: int = 0):
    import scanpy as sc
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score

    a = adata.copy()
    # Use X_integrated if present else X_uncorrected_pca
    if "X_integrated" in a.obsm:
        use_rep = "X_integrated"
    else:
        use_rep = "X_uncorrected_pca"

    # Map __hidden_label__ → label (reveal for metric)
    true_labels = a.obs.get("__hidden_label__", a.obs.get(label_key, pd.Series(["?"] * a.n_obs))).astype(str).to_numpy()

    # Leiden clustering on integrated embedding
    sc.pp.neighbors(a, use_rep=use_rep, random_state=seed, n_neighbors=15)
    sc.tl.leiden(a, resolution=1.0, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
    leiden = a.obs["leiden"].astype(int).to_numpy()

    ari = float(adjusted_rand_score(true_labels, leiden))
    nmi = float(normalized_mutual_info_score(true_labels, leiden))

    # Silhouette: use a subsample to keep it fast
    rng = np.random.default_rng(seed)
    if a.n_obs > 5000:
        idx = rng.choice(a.n_obs, size=5000, replace=False)
    else:
        idx = np.arange(a.n_obs)
    emb = a.obsm[use_rep][idx]
    lab = true_labels[idx]
    # Drop classes with fewer than 2 members
    uniq, counts = np.unique(lab, return_counts=True)
    keep_lab = set(uniq[counts >= 2])
    mask = np.array([l in keep_lab for l in lab])
    if mask.sum() > 2 and len(set(lab[mask])) >= 2:
        asw = float((silhouette_score(emb[mask], lab[mask]) + 1) / 2)  # normalize to [0,1]
    else:
        asw = float("nan")

    # Graph connectivity: per label, fraction of cells in the largest connected component.
    from scipy.sparse.csgraph import connected_components
    from scipy.sparse import csr_matrix
    conn_per_label = []
    for l in keep_lab:
        sel = np.where(true_labels == l)[0]
        if len(sel) < 3:
            continue
        # Build adjacency on the k-NN graph restricted to label cells
        knn = a.obsp["connectivities"][sel][:, sel]
        knn = csr_matrix(knn > 0, dtype=bool)
        _, labels_cc = connected_components(knn, directed=False)
        largest = max(np.bincount(labels_cc)) if len(labels_cc) else 0
        conn_per_label.append(largest / max(len(sel), 1))
    graph_conn = float(np.mean(conn_per_label)) if conn_per_label else float("nan")

    # iLISI (batch mixing) and cLISI (label separation) — simplified proxies.
    # Use per-cell fraction-of-same-batch and fraction-of-same-label in k-NN.
    batches = a.obs[batch_key].astype(str).to_numpy()
    n_sample = min(2000, a.n_obs)
    idx_s = rng.choice(a.n_obs, size=n_sample, replace=False)
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=51, metric="euclidean", n_jobs=-1)
    nn.fit(a.obsm[use_rep])
    _, knn_idx = nn.kneighbors(a.obsm[use_rep][idx_s])
    knn_idx = knn_idx[:, 1:]  # drop self
    # iLISI: number of unique batches in kNN normalized by max possible
    ilisi_values = []
    clisi_values = []
    for i, ni in enumerate(knn_idx):
        nb = batches[ni]
        ilisi_values.append(len(set(nb)) / min(len(set(batches)), 50))
        nl = true_labels[ni]
        clisi_values.append(1.0 - len(set(nl)) / 50)
    ilisi = float(np.mean(ilisi_values))
    clisi = float(np.mean(clisi_values))

    return {
        "ari": ari,
        "nmi": nmi,
        "asw": asw,
        "graph_conn": graph_conn,
        "ilisi": ilisi,
        "clisi": clisi,
    }


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    import anndata as ad

    out_dir = SHARED / "scib"
    out_dir.mkdir(parents=True, exist_ok=True)

    targets = [
        # (method, dataset, filename stem)
        ("harmony", "retina_bc89_holdout"),
        ("scanorama", "retina_bc89_holdout"),
        ("scvi", "retina_bc89_holdout"),
        ("harmony", "cheng_paad_lamp3_holdout"),
        ("scanorama", "cheng_paad_lamp3_holdout"),
        ("scvi", "cheng_paad_lamp3_holdout"),
        ("harmony", "cheng5_lamp3_holdout"),
        ("scanorama", "cheng5_lamp3_holdout"),
        ("scvi", "cheng5_lamp3_holdout"),
        ("harmony", "pancreas_synth"),
        ("scanorama", "pancreas_synth"),
        ("scvi", "pancreas_synth"),
    ]

    results = []
    for method, dataset in targets:
        h5ad = SHARED / "integrations" / method / f"{dataset}.h5ad"
        if not h5ad.exists():
            logger.warning("skip %s/%s (no H5AD)", method, dataset)
            continue
        logger.info("computing scIB for %s/%s", method, dataset)
        a = ad.read_h5ad(h5ad)
        try:
            metrics = compute_scib(a)
            out = out_dir / f"{method}_{dataset}.json"
            out.write_text(json.dumps(metrics, indent=2))
            logger.info("  %s", metrics)
            results.append({"method": method, "dataset": dataset, **metrics})
        except Exception as exc:
            logger.exception("%s/%s failed: %s", method, dataset, exc)

    df = pd.DataFrame(results)
    tsv = out_dir / "all.tsv"
    df.to_csv(tsv, sep="\t", index=False)
    logger.info("wrote %s", tsv)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
