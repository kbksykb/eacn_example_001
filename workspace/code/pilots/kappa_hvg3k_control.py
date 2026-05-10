"""Discriminating Candidate A (HVG3K preprocessing) vs Candidate B (state-regime
mechanism) for the Sade κ disconfirmation.

Run Scanorama + scVI on Cheng6_lamp3_holdout using HVG3K preprocessing (matching
the Sade pipeline). Compute κ on those integrations. Compare to the full-gene
Cheng6 κ values:

  Expected under Candidate A: κ drops significantly (Sade-like levels), because
  HVG subsetting shrinks centroid-displacement.
  Expected under Candidate B: κ stays high (Cheng-like), because the mechanism
  is cluster-regime vs state-regime, not preprocessing-dependent.
"""

from __future__ import annotations

from pathlib import Path
import time

import anndata as ad
import numpy as np
import scanpy as sc


SHARED = Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared")


def hvg_subset_cheng6(n_hvg: int = 3000) -> ad.AnnData:
    in_path = SHARED / "data" / "cheng6_lamp3_holdout.h5ad"
    a = ad.read_h5ad(in_path)
    print(f"[hvg] {in_path.name} shape={a.shape}")
    sc.pp.highly_variable_genes(a, n_top_genes=n_hvg, batch_key="batch", flavor="seurat_v3")
    a_hvg = a[:, a.var["highly_variable"]].copy()
    print(f"  hvg{n_hvg} subset shape {a_hvg.shape}")
    return a_hvg


def scanorama_on_hvg(a_hvg):
    import scanorama
    import pandas as pd
    batches = sorted(a_hvg.obs["batch"].astype(str).unique())
    print(f"[scanorama] {len(batches)} batches")
    datasets = [a_hvg[a_hvg.obs["batch"].astype(str) == b, :] for b in batches]
    data_list = [d.X.toarray() if hasattr(d.X, "toarray") else np.asarray(d.X) for d in datasets]
    gene_list = [d.var_names.tolist() for d in datasets]
    t0 = time.time()
    corrected, _ = scanorama.correct(data_list, gene_list, dimred=32)
    print(f"  scanorama done in {time.time()-t0:.0f}s; first shape {corrected[0].shape}")
    # reassemble
    out = np.zeros((a_hvg.n_obs, corrected[0].shape[1]), dtype=np.float32)
    pos = 0
    for d, c in zip(datasets, corrected):
        n = d.n_obs
        out[a_hvg.obs_names.get_indexer(d.obs_names)] = np.asarray(c, dtype=np.float32)
        pos += n
    return out


def scvi_on_hvg(a_hvg):
    import scvi
    import torch
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(a_hvg, batch_key="batch")
    m = scvi.model.SCVI(a_hvg, n_latent=30)
    print(f"[scvi] torch cuda={torch.cuda.is_available()}; training max_epochs=200")
    m.train(max_epochs=200, early_stopping=True, check_val_every_n_epoch=10)
    return m.get_latent_representation().astype(np.float32)


def per_motif_kappa_from_arrays(z_pre: np.ndarray, z_post: np.ndarray, motif_ids: np.ndarray, batch: np.ndarray) -> tuple[float, float, float]:
    kappas = []
    for cid in np.unique(motif_ids):
        mask = motif_ids == cid
        if mask.sum() < 10:
            continue
        batches_in_motif = np.unique(batch[mask])
        if len(batches_in_motif) < 2:
            continue
        dist_sum = 0.0
        n_pairs = 0
        for i, bi in enumerate(batches_in_motif):
            for bj in batches_in_motif[i+1:]:
                mi_pre = z_pre[mask & (batch == bi)].mean(axis=0) if (mask & (batch == bi)).sum() > 0 else None
                mj_pre = z_pre[mask & (batch == bj)].mean(axis=0) if (mask & (batch == bj)).sum() > 0 else None
                mi_post = z_post[mask & (batch == bi)].mean(axis=0) if (mask & (batch == bi)).sum() > 0 else None
                mj_post = z_post[mask & (batch == bj)].mean(axis=0) if (mask & (batch == bj)).sum() > 0 else None
                if mi_pre is None or mj_pre is None or mi_post is None or mj_post is None:
                    continue
                pre_sep = float(np.linalg.norm(mi_pre - mj_pre))
                post_sep = float(np.linalg.norm(mi_post - mj_post))
                dist_sum += abs(pre_sep - post_sep)
                n_pairs += 1
        if n_pairs > 0:
            kappas.append(dist_sum / n_pairs)
    if not kappas:
        return 0.0, 0.0, 0.0
    k = np.asarray(kappas, dtype=np.float32)
    return float(np.median(k)), float(np.percentile(k, 25)), float(np.percentile(k, 75))


def main() -> None:
    a_hvg = hvg_subset_cheng6(n_hvg=3000)

    # Get pre-integration PCA + motif_id_pre from the FULL-GENE Harmony h5ad
    ref = ad.read_h5ad(SHARED / "integrations" / "harmony" / "cheng6_lamp3_holdout.h5ad", backed="r")
    assert list(ref.obs_names) == list(a_hvg.obs_names), "obs_names mismatch"
    motif_ids = ref.obs["motif_id_pre"].astype(int).to_numpy()
    import pandas as pd
    batch = pd.Categorical(a_hvg.obs["batch"]).codes.astype(np.int32)

    # For pre-integration, recompute PCA on HVG3K subset for consistency
    print("[pre] recomputing PCA on HVG3K subset")
    sc.pp.pca(a_hvg, n_comps=32, random_state=0)
    z_pre = a_hvg.obsm["X_pca"].astype(np.float32)

    import json
    out_lines = []

    # Scanorama
    try:
        z_post = scanorama_on_hvg(a_hvg)
        km, kq1, kq3 = per_motif_kappa_from_arrays(z_pre, z_post, motif_ids, batch)
        line = f"cheng6_hvg3k scanorama kappa_median={km:.3f} IQR=[{kq1:.3f},{kq3:.3f}]  vs full-gene Cheng6 scanorama kappa=9.13"
        print(line); out_lines.append(line)
    except Exception as e:
        print(f"[scanorama] ERR: {e}")
        out_lines.append(f"cheng6_hvg3k scanorama ERR {e}")

    # scVI
    try:
        z_post = scvi_on_hvg(a_hvg)
        km, kq1, kq3 = per_motif_kappa_from_arrays(z_pre, z_post, motif_ids, batch)
        line = f"cheng6_hvg3k scvi kappa_median={km:.3f} IQR=[{kq1:.3f},{kq3:.3f}]  vs full-gene Cheng6 scvi kappa=5.13"
        print(line); out_lines.append(line)
    except Exception as e:
        print(f"[scvi] ERR: {e}")
        out_lines.append(f"cheng6_hvg3k scvi ERR {e}")

    out_path = Path("workspace/results/hvg3k_vs_full_kappa_control.txt")
    out_path.write_text("\n".join(out_lines) + "\n")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
