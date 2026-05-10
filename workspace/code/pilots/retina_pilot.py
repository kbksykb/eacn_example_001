"""
End-to-end REAL pilot on the Shekhar 2016 retinal bipolar atlas.

Dataset: scvi.data.retina() — 19,829 cells, 13,166 genes, 2 batches, 15 types.
Rare positive controls: BC8_9 (~1.1%), BC4 (~1.5%), BC5B (~1.8%), BC2 / BC3A / BC3B (~2-3%).

Holdout strategy (hide-and-measure): mask BC8_9 (rarest) label, run integrations,
score REAL, compare to held-out BC8_9 membership post-hoc. AUPRC is the primary
metric (per ML's RareShield validation gate).
"""

from __future__ import annotations

import argparse
import logging
import pathlib
import time
import json

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

SHARED = pathlib.Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", default=str(SHARED / "data" / "retina_scvi.h5ad"))
    p.add_argument("--rare-label", default="BC8_9")
    p.add_argument("--methods", nargs="+", default=["harmony", "scanorama", "scvi"])
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out-name", default="retina_bc89_holdout")
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")

    import anndata as ad
    import scanpy as sc
    from workspace.code.real_channels import ChannelConfig, loss_rate_at_k, score_purity_and_procrustes
    from workspace.code.real_channels import ot_wrapper
    from workspace.code.real_channels.ot_channel import OTChannelConfig
    from workspace.code.pilots.synth_pilot import (
        METHODS, _pca_project, _peak_ram_gb, _peak_vram_gb, _reset_vram,
        _git_sha, over_cluster, preprocess,
    )

    logger.info("loading retina H5AD %s", args.h5ad)
    adata = ad.read_h5ad(args.h5ad)
    adata.obs["cell_type_public"] = adata.obs["labels"].astype("category")
    adata.obs["batch"] = adata.obs["batch"].astype(str).astype("category")

    # Mask the rare label
    is_rare = (adata.obs["labels"].astype(str) == args.rare_label).to_numpy()
    adata.obs["__hidden_label__"] = adata.obs["labels"].astype(str).copy()
    adata.obs["heldout_rare"] = is_rare
    adata.obs["oracle_rare_score"] = np.where(is_rare, 2.0, 0.0).astype(np.float32)
    new = adata.obs["cell_type_public"].astype(object).copy()
    new[is_rare] = np.nan
    adata.obs["cell_type_public"] = new.astype("category")
    logger.info("rare cells (%s) masked: %d", args.rare_label, int(is_rare.sum()))

    # Preprocess
    adata_pre = preprocess(adata, batch_key="batch", seed=args.seed)
    candidates = over_cluster(adata_pre, seed=args.seed, resolution=2.0)
    adata_pre.obs["motif_id_pre"] = candidates
    logger.info("Leiden motifs: %d", len(np.unique(candidates)))

    # Ground-truth rare-truth motif ids
    hidden = adata.obs["__hidden_label__"].astype(str).to_numpy()
    rare_ids: set[int] = set()
    is_rare = np.asarray([x == args.rare_label for x in hidden])
    for cid in np.unique(candidates):
        m = candidates == cid
        if m.sum() == 0:
            continue
        if is_rare[m].sum() >= 3 and is_rare[m].sum() / m.sum() > 0.5:
            rare_ids.add(int(cid))
    logger.info("rare-truth motif ids (overlap > 50%%, ≥3 cells): %s", sorted(rare_ids))

    run_info = {
        "random_state": args.seed,
        "commit_sha": _git_sha(),
        "run_date": time.strftime("%Y-%m-%d"),
    }
    adata_pre.uns["reproducibility"] = run_info

    # Write pre
    (SHARED / "data").mkdir(parents=True, exist_ok=True)
    adata_pre.write_h5ad(SHARED / "data" / f"{args.out_name}.h5ad")

    batch_int = pd.Categorical(adata_pre.obs["batch"]).codes.astype(np.int32)
    manifest_rows = []
    summary_rows = []
    det_dir = SHARED / "detections" / "real"
    det_dir.mkdir(parents=True, exist_ok=True)

    from sklearn.metrics import average_precision_score

    for method in args.methods:
        if method not in METHODS:
            logger.warning("%s not in METHODS", method)
            continue
        logger.info("=== %s ===", method)
        _reset_vram()
        t0 = time.perf_counter()
        try:
            emb_int, method_meta = METHODS[method](adata_pre, batch_key="batch", seed=args.seed)
        except Exception:
            logger.exception("%s failed", method)
            continue
        elapsed = time.perf_counter() - t0
        vram = _peak_vram_gb()
        ram = _peak_ram_gb()

        # Save post
        out = adata_pre.copy()
        out.obsm["X_integrated"] = emb_int
        out.uns["method"] = method
        out.uns["method_version"] = method_meta.get("method_version", "unknown")
        out.uns["method_config"] = {"seed": args.seed, "d_projection": 32}
        out.uns["runtime"] = {
            "wall_clock_seconds": round(elapsed, 3),
            "peak_vram_gb": round(vram, 3),
            "peak_ram_gb": round(ram, 3),
            "n_cells_used": int(adata_pre.n_obs),
        }
        (SHARED / "integrations" / method).mkdir(parents=True, exist_ok=True)
        out.write_h5ad(SHARED / "integrations" / method / f"{args.out_name}.h5ad")

        # REAL channels
        d = emb_int.shape[1]
        pre_trunc = adata_pre.obsm["X_uncorrected_pca"][:, :d].astype(np.float32)
        cfg = ChannelConfig(k=15, n_bootstrap=20, rare_abundance_threshold=0.05, seed=args.seed)
        reports = score_purity_and_procrustes(
            emb_pre=pre_trunc,
            emb_post=emb_int.astype(np.float32),
            candidate_labels=candidates,
            batch_labels=batch_int,
            cfg=cfg,
        )

        # OT channel
        ot_cfg = OTChannelConfig(n_permutations=100, device="cpu", k_neighbors=15, seed=args.seed)
        ot_res = ot_wrapper.score_two_sided(
            emb_pre=pre_trunc,
            emb_post=emb_int.astype(np.float32),
            candidate_labels=candidates,
            batch_labels=batch_int,
            cfg=ot_cfg,
        )
        ot_by_id = {
            int(cid): (stat, p, s) for cid, stat, p, s in zip(
                ot_res["candidate_ids"], ot_res["stat"],
                ot_res["p_values"], ot_res["signed_stat_mean"]
            )
        }

        rows = []
        y_true = []
        y_score = []
        for r in reports:
            stat, p, signed = ot_by_id.get(r.candidate_id, (np.nan, np.nan, np.nan))
            gt = r.candidate_id in rare_ids
            rows.append({
                "method": method,
                "dataset": args.out_name,
                "candidate_id": r.candidate_id,
                "abundance": r.abundance,
                "channel_mknn_purity_pre": r.purity_pre,
                "channel_mknn_purity_post": r.purity_post,
                "channel_proc_displacement": r.procrustes_displacement,
                "channel_boot_stable": r.bootstrap_fraction_stable,
                "channel_ot_stat": stat,
                "channel_ot_p_value": p,
                "channel_ot_signed": signed,
                "loss_probability": r.loss_probability,
                "ground_truth_label": gt,
                "seed": args.seed,
                "n_cells": int(adata_pre.n_obs),
            })
            y_true.append(bool(gt))
            y_score.append(r.loss_probability)

        df = pd.DataFrame(rows)
        df.to_parquet(det_dir / f"{method}_{args.out_name}.parquet")

        lr1 = loss_rate_at_k(reports, k=1)
        lr3 = loss_rate_at_k(reports, k=3)
        y_true_arr = np.asarray(y_true)
        y_score_arr = np.asarray(y_score)
        auprc = float(average_precision_score(y_true_arr, y_score_arr)) if y_true_arr.any() else float("nan")

        logger.info(
            "%s: parquet written, LossRate@1=%.3f, @3=%.3f, AUPRC=%.3f, wall=%.1fs",
            method, lr1, lr3, auprc, elapsed,
        )

        # Min OT p on rare-truth and non-rare-truth
        min_p_rare = float(df.loc[df["ground_truth_label"], "channel_ot_p_value"].min()) if df["ground_truth_label"].any() else float("nan")
        min_p_nonrare = float(df.loc[~df["ground_truth_label"], "channel_ot_p_value"].min()) if (~df["ground_truth_label"]).any() else float("nan")
        summary_rows.append({
            "method": method,
            "LossRate@1": lr1,
            "LossRate@3": lr3,
            "AUPRC": auprc,
            "min_ot_p_rare": min_p_rare,
            "min_ot_p_nonrare": min_p_nonrare,
            "n_candidates": len(reports),
            "n_rare_truth_candidates": len(rare_ids),
        })

        manifest_rows.append({
            "method": method,
            "dataset": args.out_name,
            "wall_clock_seconds": round(elapsed, 3),
            "peak_vram_gb": round(vram, 3),
            "peak_ram_gb": round(ram, 3),
            "n_cells_used": int(adata_pre.n_obs),
            "package_version": method_meta.get("method_version", "unknown"),
            "commit_sha": run_info["commit_sha"],
            "run_date": run_info["run_date"],
        })

    manifest_path = SHARED / "integrations" / "manifest.csv"
    mdf = pd.DataFrame(manifest_rows)
    if manifest_path.exists():
        existing = pd.read_csv(manifest_path)
        existing = existing[~((existing["method"].isin(mdf["method"])) & (existing["dataset"] == args.out_name))]
        mdf = pd.concat([existing, mdf], ignore_index=True)
    mdf.to_csv(manifest_path, index=False)

    print("\n=== RETINA PILOT SUMMARY (BC8_9 holdout) ===")
    sdf = pd.DataFrame(summary_rows)
    print(sdf.to_string(index=False))
    results_dir = pathlib.Path(__file__).resolve().parents[3] / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    sdf.to_csv(results_dir / f"{args.out_name}_summary.csv", index=False)


if __name__ == "__main__":
    main()
