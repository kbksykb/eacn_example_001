"""BBKNN blind prospective run on Cheng PAAD + Cheng6 × LAMP3+ mregDC.

Pre-registered prediction (workspace/results/bbknn_prospective_prediction.md):
  κ_med ≤ 2.5, Regime A (scale-monotonic).
  Cheng-PAAD p_bh ≥ 0.5; Cheng6 p_bh < 0.01.

This runner is BLIND — writes integrations + detections + κ but does NOT
read/print them. Results inspected separately by bbknn_prospective_eval.py.

Uses Harmony's pipeline output as source of motif_id_pre + __hidden_label__
(BBKNN operates on neighbors-graph and emits X_integrated via PCA-align-style
averaging; we reuse the same pre-integration PCA for consistency).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

SHARED = Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared")


def bbknn_integrate(a: ad.AnnData) -> np.ndarray:
    """Run BBKNN and return 32-dim 'integrated' embedding.

    BBKNN returns a batch-balanced kNN graph (anndata.uns['neighbors']), not a
    direct embedding. To get a point embedding we compute diffmap on the BBKNN
    graph and take 32 components. This is the standard BBKNN → REAL adapter.
    """
    import bbknn
    # BBKNN expects X_pca with standard PCA components
    sc.pp.pca(a, n_comps=50, random_state=0)
    t0 = time.time()
    bbknn.bbknn(a, batch_key="batch", n_pcs=50, neighbors_within_batch=3, trim=50)
    print(f"  bbknn graph done in {time.time()-t0:.0f}s")
    # Compute a 32-dim representation via PCA on the BBKNN-connectivity-weighted X_pca
    # Simpler: use the rotated 50-dim PCA as X_integrated (BBKNN doesn't modify features,
    # only the graph). For κ + REAL scoring this captures the geometric change because
    # BBKNN effectively re-weights neighbors rather than transforming coordinates.
    # Proper alternative: diffmap(a, n_comps=32) on the BBKNN graph.
    sc.tl.diffmap(a, n_comps=32)
    Z = a.obsm["X_diffmap"].astype(np.float32)
    print(f"  diffmap embedding shape {Z.shape}")
    return Z


def run_on_holdout(holdout_name: str) -> Path:
    input_ref = SHARED / "integrations" / "harmony" / f"{holdout_name}.h5ad"
    out_path = SHARED / "integrations" / "bbknn" / f"{holdout_name}.h5ad"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"[{time.strftime('%H:%M:%S')}] loading source from harmony {input_ref.name}")
    ref = ad.read_h5ad(input_ref)
    # Build a NEW working AnnData from the holdout pre-integration data
    input_holdout = SHARED / "data" / f"{holdout_name}.h5ad"
    a = ad.read_h5ad(input_holdout)
    # Copy pipeline-critical obs columns from ref to preserve motif_id_pre + __hidden_label__
    for key in ["motif_id_pre", "__hidden_label__", "heldout_rare", "cell_type_public", "leiden", "oracle_rare_score", "MajorCluster", "labels"]:
        if key in ref.obs and key not in a.obs:
            a.obs[key] = ref.obs[key].values
    # Ensure pre-integration PCA is present
    if "X_uncorrected_pca" in ref.obsm:
        a.obsm["X_uncorrected_pca"] = ref.obsm["X_uncorrected_pca"].astype(np.float32)
    print(f"  shape {a.shape}; running BBKNN")

    Z_int = bbknn_integrate(a)
    a.obsm["X_integrated"] = Z_int

    a.write_h5ad(out_path)
    print(f"[{time.strftime('%H:%M:%S')}] wrote {out_path}")
    return out_path


def score_real_and_kappa(h5ad_path: Path, method: str, dataset: str, n: int, pi: float) -> None:
    """Run REAL OT scoring (B=500) + κ computation. BLIND — print nothing about p/κ values."""
    sys.path.insert(0, "/mnt/d-1274477442621830-m5ObBqn4/eacn3-main/examples/computational_biology")
    from workspace.code.real_channels import ot_wrapper
    from workspace.code.real_channels.ot_channel import OTChannelConfig
    from scipy.stats import false_discovery_control

    a = ad.read_h5ad(h5ad_path)
    z_pre = a.obsm["X_uncorrected_pca"].astype(np.float32)
    z_post = a.obsm["X_integrated"].astype(np.float32)
    d = z_post.shape[1]
    pre_trunc = z_pre[:, :d]
    motif_ids = a.obs["motif_id_pre"].astype(int).to_numpy()
    batch_int = pd.Categorical(a.obs["batch"]).codes.astype(np.int32)

    print(f"[{time.strftime('%H:%M:%S')}] OT B=500 on {method}/{dataset}")
    cfg = OTChannelConfig(n_permutations=500, device="cpu", k_neighbors=15, seed=0)
    r = ot_wrapper.score_two_sided(pre_trunc, z_post, motif_ids, batch_int, cfg)
    raw_p = np.asarray(r["p_values"])
    bh_p = false_discovery_control(raw_p, method="bh")

    hidden = a.obs["__hidden_label__"].astype(str).to_numpy()
    is_lamp3 = np.array([("LAMP3" in x) or ("cDC3" in x) or x == "M05_cDC3_LAMP3" for x in hidden])
    rare_ids = set()
    for cid in r["candidate_ids"]:
        mask = motif_ids == int(cid)
        if mask.sum() == 0:
            continue
        if is_lamp3[mask].sum() >= 3 and is_lamp3[mask].sum() / mask.sum() > 0.5:
            rare_ids.add(int(cid))

    rows = []
    for cid, stat, p, bp, signed in zip(r["candidate_ids"], r["stat"], raw_p, bh_p, r["signed_stat_mean"]):
        cid = int(cid)
        is_rare = cid in rare_ids
        rows.append(dict(
            method=method, dataset=dataset, candidate_id=cid,
            abundance=float((motif_ids == cid).mean()),
            channel_mknn_purity_pre=np.nan, channel_mknn_purity_post=np.nan,
            channel_proc_displacement=np.nan, channel_boot_stable=np.nan,
            channel_ot_stat=float(stat), channel_ot_p_value=float(p), channel_ot_signed=float(signed),
            loss_probability=float(np.clip(1.0 - bp, 0.0, 1.0)),
            ground_truth_label=is_rare, seed=0, n_cells=int(a.n_obs),
        ))
    df = pd.DataFrame(rows)
    det_path = SHARED / "detections" / "real" / f"{method}_{dataset}.parquet"
    det_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(det_path, index=False)
    print(f"  wrote {det_path}  rows={len(df)} (values NOT inspected)")

    # κ (per_motif_kappa inline so we don't have to import from kappa_regression)
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
    k_out = SHARED / "detections" / "real" / f"{method}_{dataset}_kappa.npy"
    np.save(k_out, np.asarray(kappas, dtype=np.float32))
    print(f"  wrote {k_out} (kappa values NOT printed; inspect via bbknn_prospective_eval.py)")


def main() -> None:
    for name, n, pi in [("cheng_paad_lamp3_holdout", 2853, 0.0107),
                         ("cheng6_lamp3_holdout", 49271, 0.0115)]:
        out = run_on_holdout(name)
        score_real_and_kappa(out, "bbknn", name, n, pi)
    print("BLIND run complete. Inspect via bbknn_prospective_eval.py.")


if __name__ == "__main__":
    main()
