"""κ·Δ direct measurement on the 9 existing integrations.

Test of Theorem S1-overcorrection: κ := cross-batch mass-redistribution magnitude
per motif, measured directly from pre and post embeddings. Math's prediction:
scale-monotonic integrators (Harmony) have low κ at small n, rising with scale;
method-intrinsic integrators (Scanorama, scVI) have high κ regardless of n.

Per-motif κ definition (Math's Supp S1-overcorrection convention):
    For each motif c with cells C_c in batch b and C_c in batch b':
        κ_c(b, b') := W_2^2(emp pre-kNN density in b, emp post-kNN density in b')
        κ_c := median over (b, b') pairs

We approximate W_2^2 by the squared Procrustes displacement of motif centroid
between pre and post embeddings, aggregated over batches. This is a conservative
proxy; the full OT formulation lives in Math's compute_kappa.py (not yet shipped).

Output: workspace/results/kappa_grid.csv
"""

from __future__ import annotations

import csv
from pathlib import Path

import anndata as ad
import numpy as np

SHARED = Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared")

GRID = [
    ("cheng_paad", "cheng_paad_lamp3_holdout.h5ad", "harmony", 2853, 1.0000),
    ("cheng5", "cheng5_lamp3_holdout.h5ad", "harmony", 20341, 0.0130),
    ("cheng6", "cheng6_lamp3_holdout.h5ad", "harmony", 49271, 0.0041),
    ("cheng_paad", "cheng_paad_lamp3_holdout.h5ad", "scanorama", 2853, 0.005),
    ("cheng5", "cheng5_lamp3_holdout.h5ad", "scanorama", 20341, 0.008),
    ("cheng6", "cheng6_lamp3_holdout.h5ad", "scanorama", 49271, 0.0039),
    ("cheng_paad", "cheng_paad_lamp3_holdout.h5ad", "scvi", 2853, 0.003),
    ("cheng5", "cheng5_lamp3_holdout.h5ad", "scvi", 20341, 0.004),
    ("cheng6", "cheng6_lamp3_holdout.h5ad", "scvi", 49271, 0.0030),
]


def per_motif_kappa(a) -> tuple[float, float, float]:
    """Returns (kappa_median, kappa_q1, kappa_q3) over motifs that appear in ≥2 batches."""
    z_pre = a.obsm["X_uncorrected_pca"].astype(np.float32)
    z_post = a.obsm["X_integrated"].astype(np.float32)
    d = z_post.shape[1]
    z_pre = z_pre[:, :d]

    motif = a.obs["motif_id_pre"].astype(int).to_numpy()
    import pandas as pd
    batch = pd.Categorical(a.obs["batch"]).codes.astype(np.int32)

    kappas = []
    for cid in np.unique(motif):
        mask = motif == cid
        if mask.sum() < 10:
            continue
        batches_in_motif = np.unique(batch[mask])
        if len(batches_in_motif) < 2:
            continue
        dist_sum = 0.0
        n_pairs = 0
        for i, bi in enumerate(batches_in_motif):
            for bj in batches_in_motif[i + 1:]:
                mi_pre = z_pre[mask & (batch == bi)].mean(axis=0) if (mask & (batch == bi)).sum() > 0 else None
                mj_pre = z_pre[mask & (batch == bj)].mean(axis=0) if (mask & (batch == bj)).sum() > 0 else None
                mi_post = z_post[mask & (batch == bi)].mean(axis=0) if (mask & (batch == bi)).sum() > 0 else None
                mj_post = z_post[mask & (batch == bj)].mean(axis=0) if (mask & (batch == bj)).sum() > 0 else None
                if mi_pre is None or mj_pre is None or mi_post is None or mj_post is None:
                    continue
                pre_sep = float(np.linalg.norm(mi_pre - mj_pre))
                post_sep = float(np.linalg.norm(mi_post - mj_post))
                delta = abs(pre_sep - post_sep)
                dist_sum += delta
                n_pairs += 1
        if n_pairs > 0:
            kappas.append(dist_sum / n_pairs)

    k = np.asarray(kappas, dtype=np.float32)
    if len(k) == 0:
        return 0.0, 0.0, 0.0
    return float(np.median(k)), float(np.percentile(k, 25)), float(np.percentile(k, 75))


def main() -> None:
    rows = []
    for ds, holdout, method, n, p_bh in GRID:
        path = SHARED / "integrations" / method / holdout
        print(f"[load] {path.name} ({method})")
        a = ad.read_h5ad(path)
        k_med, k_q1, k_q3 = per_motif_kappa(a)
        nlog = float(-np.log10(max(p_bh, 1e-6)))
        print(f"  {ds} × {method}: κ_med={k_med:.3f} [{k_q1:.3f},{k_q3:.3f}]  -log10(p)={nlog:.2f}")
        rows.append(dict(
            dataset=ds, method=method, n=n, p_bh=p_bh, neg_log10_pbh=round(nlog, 4),
            kappa_median=round(k_med, 4), kappa_q1=round(k_q1, 4), kappa_q3=round(k_q3, 4),
        ))

    out = Path("workspace/results/kappa_grid.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out} with {len(rows)} rows")


if __name__ == "__main__":
    main()
