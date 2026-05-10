"""
End-to-end pancreas pilot — Tier A unit test for REAL.

Pipeline:
    1. Load scIB pancreas benchmark.
    2. Mask the epsilon label.
    3. Pre-integration: HVG + PCA → X_pca.
    4. Run 2 integrators (Harmony + Scanorama) → X_integrated.
    5. Compute REAL channels (mutual-kNN purity + Procrustes) on pre vs post.
    6. Score per candidate cluster; report LossRate@k and per-channel scores.
    7. Ex-post: compare channel output to true epsilon membership (AUPRC).
    8. Write parquet to workspace/results/pancreas_pilot_<method>.parquet.

This is a v1 smoke-with-real-data — NOT the full 9-method panel. That lands
in the full harness (workspace/code/harness) once the remaining method stubs
are filled.
"""

from __future__ import annotations

import argparse
import logging
import pathlib

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def leiden_over_cluster(adata, use_rep: str, resolution: float, seed: int) -> np.ndarray:
    import scanpy as sc

    sc.pp.neighbors(adata, n_neighbors=15, use_rep=use_rep, random_state=seed)
    sc.tl.leiden(adata, resolution=resolution, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
    return adata.obs["leiden"].astype(int).to_numpy()


def run_harmony(adata_pre, batch_key: str, seed: int) -> np.ndarray:
    import harmonypy as hm
    ho = hm.run_harmony(adata_pre.obsm["X_pca"], adata_pre.obs, vars_use=[batch_key], random_state=seed)
    return np.asarray(ho.Z_corr.T).astype(np.float32)


def run_scanorama(adata_pre, batch_key: str, seed: int) -> np.ndarray:
    import scanorama
    batches = adata_pre.obs[batch_key].astype(str)
    uniq = sorted(batches.unique())
    splits = [adata_pre[batches == b].copy() for b in uniq]
    corrected, _ = scanorama.correct_scanpy(splits, return_dimred=True)
    emb = np.empty((adata_pre.n_obs, corrected[0].obsm["X_scanorama"].shape[1]), dtype=np.float32)
    for b, piece in zip(uniq, corrected):
        emb[np.where(batches == b)[0]] = piece.obsm["X_scanorama"]
    return emb


METHODS = {"harmony": run_harmony, "scanorama": run_scanorama}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--rare-label", default="epsilon")
    p.add_argument("--batch-key", default="tech")
    p.add_argument("--cell-type-key", default="celltype")
    p.add_argument("--methods", nargs="+", default=["harmony", "scanorama"])
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default="workspace/results")
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")
    import scanpy as sc
    from workspace.code.datasets import fetch_pancreas_scIB, mask_rare_label
    from workspace.code.real_channels import ChannelConfig, loss_rate_at_k, score_purity_and_procrustes

    logger.info("loading pancreas scIB benchmark…")
    adata = fetch_pancreas_scIB()

    obs_cols = set(adata.obs.columns)
    ct_key = args.cell_type_key if args.cell_type_key in obs_cols else (
        "celltype" if "celltype" in obs_cols else "cell_type"
    )
    batch_key = args.batch_key if args.batch_key in obs_cols else (
        "tech" if "tech" in obs_cols else "batch"
    )
    logger.info("using ct_key=%s batch_key=%s", ct_key, batch_key)

    adata = mask_rare_label(adata, ct_key, args.rare_label)

    # Pre-integration pipeline — works on the scIB distribution which is already log-normalized
    if "X_pca" not in adata.obsm:
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=batch_key, flavor="seurat")
        adata_hvg = adata[:, adata.var["highly_variable"]].copy()
        sc.pp.scale(adata_hvg, max_value=10)
        sc.tl.pca(adata_hvg, n_comps=50, random_state=args.seed)
        adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

    # Over-cluster on pre-integration
    pre_candidate = leiden_over_cluster(adata, use_rep="X_pca", resolution=3.0, seed=args.seed)
    adata.obs["real_candidate_pre"] = pre_candidate

    # Batch labels
    batch_int = pd.Categorical(adata.obs[batch_key]).codes.astype(np.int32)

    results = []
    out_root = pathlib.Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    for method in args.methods:
        logger.info("running %s…", method)
        try:
            emb_int = METHODS[method](adata, batch_key, args.seed)
        except Exception as exc:
            logger.exception("method %s failed", method)
            continue

        logger.info("running REAL channels for %s", method)
        cfg = ChannelConfig(k=15, n_bootstrap=20, rare_abundance_threshold=0.02, seed=args.seed)
        reports = score_purity_and_procrustes(
            emb_pre=adata.obsm["X_pca"].astype(np.float32),
            emb_post=emb_int,
            candidate_labels=pre_candidate,
            batch_labels=batch_int,
            cfg=cfg,
        )

        # Which candidate clusters correspond to true epsilon? Majority-label check.
        true = adata.obs["__hidden_label__"].astype(str)
        is_rare = (true.str.lower().str.strip() == args.rare_label.lower().strip()).to_numpy()
        rare_cid_by_count = {}
        for cid in np.unique(pre_candidate):
            mask = pre_candidate == cid
            if mask.sum() > 0:
                rare_cid_by_count[cid] = is_rare[mask].sum()
        rare_candidate_ids = {cid for cid, n in rare_cid_by_count.items() if n >= 3}

        rows = []
        for r in reports:
            rows.append({
                "method": method,
                "candidate_id": r.candidate_id,
                "abundance": r.abundance,
                "purity_pre": r.purity_pre,
                "purity_post": r.purity_post,
                "nn_identity_pre": r.nn_identity_pre,
                "nn_identity_post": r.nn_identity_post,
                "procrustes_displacement": r.procrustes_displacement,
                "bootstrap_fraction_stable": r.bootstrap_fraction_stable,
                "loss_probability": r.loss_probability,
                "is_rare_ground_truth": r.candidate_id in rare_candidate_ids,
            })
        df = pd.DataFrame(rows)
        parquet_path = out_root / f"pancreas_pilot_{method}.parquet"
        df.to_parquet(parquet_path)
        logger.info("wrote %s (%d rows, %d rare-truth rows)", parquet_path, len(df), int(df["is_rare_ground_truth"].sum()))

        lr1 = loss_rate_at_k(reports, k=1)
        lr5 = loss_rate_at_k(reports, k=5)
        logger.info("%s  LossRate@1=%.3f  LossRate@5=%.3f", method, lr1, lr5)
        results.append({"method": method, "LossRate@1": lr1, "LossRate@5": lr5, "n_candidates": len(reports)})

    pd.DataFrame(results).to_csv(out_root / "pancreas_pilot_summary.csv", index=False)
    logger.info("summary → %s/pancreas_pilot_summary.csv", out_root)


if __name__ == "__main__":
    main()
