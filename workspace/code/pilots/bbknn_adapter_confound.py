"""Alternative BBKNN adapter: X_pca-unchanged κ test.

Under Math's intuition that BBKNN preserves kNN → low κ, the relevant κ should
be measured on the embedding BBKNN actually touches. BBKNN only modifies the
neighbor graph, not X_pca. So the κ-on-X_pca should be ≈ 0 (identity).

This is a trivial check that exposes the adapter-choice confound: if κ on X_pca
≈ 0, then our diffmap κ ≈ 9 is entirely adapter-induced, not BBKNN-induced.

The κ-on-X_integrated framework applies to integrators that produce a point-
embedding. BBKNN-graph-only has no such native embedding; users must choose an
adapter, and that choice itself is a κ-affecting variable.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

SHARED = Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared")


def per_motif_kappa(z_pre, z_post, motif_ids, batch_int):
    kappas = []
    for cid in np.unique(motif_ids):
        mask = motif_ids == cid
        if mask.sum() < 10:
            continue
        batches_in_motif = np.unique(batch_int[mask])
        if len(batches_in_motif) < 2:
            continue
        dist_sum = 0.0
        n_pairs = 0
        for i, bi in enumerate(batches_in_motif):
            for bj in batches_in_motif[i + 1:]:
                mi_pre = z_pre[mask & (batch_int == bi)].mean(axis=0) if (mask & (batch_int == bi)).sum() > 0 else None
                mj_pre = z_pre[mask & (batch_int == bj)].mean(axis=0) if (mask & (batch_int == bj)).sum() > 0 else None
                mi_post = z_post[mask & (batch_int == bi)].mean(axis=0) if (mask & (batch_int == bi)).sum() > 0 else None
                mj_post = z_post[mask & (batch_int == bj)].mean(axis=0) if (mask & (batch_int == bj)).sum() > 0 else None
                if mi_pre is None or mj_pre is None or mi_post is None or mj_post is None:
                    continue
                pre_sep = float(np.linalg.norm(mi_pre - mj_pre))
                post_sep = float(np.linalg.norm(mi_post - mj_post))
                dist_sum += abs(pre_sep - post_sep)
                n_pairs += 1
        if n_pairs > 0:
            kappas.append(dist_sum / n_pairs)
    k = np.asarray(kappas, dtype=np.float32)
    if len(k) == 0:
        return 0.0, 0.0, 0.0
    return float(np.median(k)), float(np.percentile(k, 25)), float(np.percentile(k, 75))


def main() -> None:
    rows = []
    for name in ["cheng_paad_lamp3_holdout", "cheng6_lamp3_holdout"]:
        a = ad.read_h5ad(SHARED / "integrations" / "bbknn" / f"{name}.h5ad")
        z_pre = a.obsm["X_uncorrected_pca"].astype(np.float32)
        d = z_pre.shape[1]
        # X_pca "passthrough" — BBKNN doesn't touch PCA
        z_post_pca = a.obsm["X_pca"][:, :d].astype(np.float32)
        z_post_diff = a.obsm["X_integrated"].astype(np.float32)  # what we used in blind test
        d_diff = z_post_diff.shape[1]
        motif_ids = a.obs["motif_id_pre"].astype(int).to_numpy()
        batch_int = pd.Categorical(a.obs["batch"]).codes.astype(np.int32)

        km_pca, q1_pca, q3_pca = per_motif_kappa(z_pre, z_post_pca, motif_ids, batch_int)
        km_diff, q1_diff, q3_diff = per_motif_kappa(z_pre[:, :d_diff], z_post_diff, motif_ids, batch_int)

        print(f"{name}: κ_pca(passthrough)={km_pca:.3f} IQR=[{q1_pca:.3f},{q3_pca:.3f}] | κ_diffmap={km_diff:.3f} IQR=[{q1_diff:.3f},{q3_diff:.3f}]")
        rows.append((name, km_pca, q1_pca, q3_pca, km_diff, q1_diff, q3_diff))

    out = Path("workspace/results/bbknn_adapter_confound.txt")
    with out.open("w") as f:
        f.write("BBKNN adapter confound test:\n")
        f.write("=================================================\n")
        for name, kp, qp1, qp3, kd, qd1, qd3 in rows:
            f.write(f"\n{name}\n")
            f.write(f"  κ with X_pca (passthrough):    median={kp:.3f} IQR=[{qp1:.3f},{qp3:.3f}]\n")
            f.write(f"  κ with X_diffmap (our choice): median={kd:.3f} IQR=[{qd1:.3f},{qd3:.3f}]\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
