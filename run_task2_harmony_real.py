#!/usr/bin/env python3
"""
Task 2 — REAL evaluation pipeline (synth_pilot spec).

Independent end-to-end flow: preprocess (1000 HVG) → Leiden motifs → Harmony (32-d)
→ REAL three-channel detection. Does NOT read task 1 output.

Outputs live under <out-root>/task2/ (separate from task 1).

Usage (from repo root):
    python run_task2_harmony_real.py --input workspace/Galaxy1.h5ad
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd

from pipeline_common import (
    add_common_input_args,
    dataset_name_from_input,
    git_sha_short,
    resolve_input_path,
    resolve_out_root,
    validate_batch_key,
)
from workspace.code.pilots.synth_pilot import (
    METHODS,
    over_cluster,
    preprocess,
    run_real_channels,
)
logger = logging.getLogger(__name__)


def _check_task2_dependencies() -> None:
    missing = []
    try:
        import torch  # noqa: F401
    except ImportError:
        missing.append("torch (required for REAL OT channel)")
    try:
        import igraph  # noqa: F401
    except ImportError:
        missing.append("igraph + leidenalg (required for Leiden clustering)")
    if missing:
        raise ImportError(
            "Task 2 missing dependencies: "
            + "; ".join(missing)
            + ". Install with: pip install torch igraph leidenalg pyarrow"
        )


def main() -> int:
    p = argparse.ArgumentParser(
        description="Task 2: synth_pilot REAL pipeline (preprocess + Harmony + REAL)."
    )
    add_common_input_args(p)
    p.add_argument(
        "--leiden-resolution",
        type=float,
        default=3.0,
        help="Leiden resolution for motif candidates (default: 3.0).",
    )
    p.add_argument(
        "--cell-type-key",
        default=None,
        help="Optional obs column for ex-post ground_truth_label in parquet "
        "(does not affect REAL scoring).",
    )
    p.add_argument(
        "--rare-label",
        default=None,
        help="If --cell-type-key is set, cells matching this label are ground-truth rare motifs.",
    )
    args = p.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    _check_task2_dependencies()

    input_path = resolve_input_path(args.input)
    name = dataset_name_from_input(input_path, args.name)
    out_root = resolve_out_root(args.out_root)
    task_root = out_root / "task2"

    logger.info("Task 2 — Harmony + REAL (synth_pilot)")
    logger.info("  input: %s", input_path)
    logger.info("  name:  %s", name)
    logger.info("  out:   %s", task_root)

    adata_raw = ad.read_h5ad(input_path)
    validate_batch_key(adata_raw, args.batch_key)

    t_all = time.perf_counter()

    logger.info("Preprocessing (1000 HVG + PCA 50)...")
    adata_pre = preprocess(adata_raw, batch_key=args.batch_key, seed=args.seed)

    logger.info("Leiden over-clustering (resolution=%s)...", args.leiden_resolution)
    candidates = over_cluster(
        adata_pre, seed=args.seed, resolution=args.leiden_resolution
    )
    adata_pre.obs["motif_id_pre"] = candidates
    logger.info("Motif candidates: %d clusters", len(np.unique(candidates)))

    pre_path = task_root / "data" / f"{name}_pre.h5ad"
    pre_path.parent.mkdir(parents=True, exist_ok=True)
    adata_pre.write_h5ad(pre_path)

    if "harmony" not in METHODS:
        logger.error("harmony not registered in synth_pilot.METHODS")
        return 1

    logger.info("Running Harmony (REAL-line spec, 32-d projection)...")
    t0 = time.perf_counter()
    emb_int, method_meta = METHODS["harmony"](
        adata_pre, batch_key=args.batch_key, seed=args.seed
    )
    harmony_s = time.perf_counter() - t0

    post_path = task_root / "integrations" / "harmony" / f"{name}.h5ad"
    post_path.parent.mkdir(parents=True, exist_ok=True)
    post = adata_pre.copy()
    post.obsm["X_integrated"] = emb_int
    post.uns["method"] = "harmony"
    post.uns["method_version"] = method_meta.get("method_version", "unknown")
    post.uns["pipeline"] = "task2_harmony_real"
    post.write_h5ad(post_path)

    logger.info("Running REAL channels...")
    t0 = time.perf_counter()
    reports, lr1, lr3 = run_real_channels(adata_pre, emb_int, candidates)
    real_s = time.perf_counter() - t0

    rare_ids: set[int] = set()
    if args.cell_type_key and args.rare_label:
        if args.cell_type_key not in adata_raw.obs.columns:
            logger.warning(
                "cell-type-key %r not found; skipping ground_truth_label",
                args.cell_type_key,
            )
        else:
            labels = adata_raw.obs[args.cell_type_key].astype(str).to_numpy()
            is_rare = np.asarray(
                [x == args.rare_label for x in labels]
            )
            for cid in np.unique(candidates):
                m = candidates == cid
                if (
                    m.sum() >= 3
                    and is_rare[m].sum() >= 3
                    and is_rare[m].sum() / m.sum() > 0.5
                ):
                    rare_ids.add(int(cid))

    rows = []
    for r in reports:
        rows.append(
            {
                "method": "harmony",
                "dataset": name,
                "candidate_id": r.candidate_id,
                "abundance": r.abundance,
                "channel_mknn_purity_pre": r.purity_pre,
                "channel_mknn_purity_post": r.purity_post,
                "channel_mknn_p_value": getattr(r, "mknn_p_value", float("nan")),
                "channel_proc_displacement": r.procrustes_displacement,
                "channel_proc_p_value": getattr(r, "proc_p_value", float("nan")),
                "channel_boot_stable": r.bootstrap_fraction_stable,
                "channel_boot_p_value": getattr(r, "boot_p_value", float("nan")),
                "channel_ot_stat": getattr(r, "ot_stat", float("nan")),
                "channel_ot_p_value": getattr(r, "ot_p_value", float("nan")),
                "channel_ot_signed": getattr(r, "ot_signed_mean", float("nan")),
                "fisher_chi2": getattr(r, "fisher_chi2", float("nan")),
                "fisher_p_value": getattr(r, "fisher_p_value", float("nan")),
                "loss_probability": r.loss_probability,
                "ground_truth_label": r.candidate_id in rare_ids,
                "seed": args.seed,
                "n_cells": int(adata_pre.n_obs),
            }
        )

    parquet_path = task_root / "detections" / "real" / f"harmony_{name}.parquet"
    parquet_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_parquet(parquet_path)

    manifest = {
        "task": "task2_harmony_real",
        "pipeline": "workspace/code/pilots/synth_pilot.py (preprocess + run_real_channels)",
        "spec": "1000 HVG, PCA 50, Harmony on X_pca → 32-d, Leiden motifs, REAL 3-channel",
        "input_h5ad": str(input_path),
        "dataset_name": name,
        "batch_key": args.batch_key,
        "seed": args.seed,
        "leiden_resolution": args.leiden_resolution,
        "n_motif_candidates": int(len(np.unique(candidates))),
        "LossRate@1": lr1,
        "LossRate@3": lr3,
        "harmony_wallclock_s": harmony_s,
        "real_wallclock_s": real_s,
        "total_wallclock_s": time.perf_counter() - t_all,
        "git_commit": git_sha_short(),
        "outputs": {
            "pre_h5ad": str(pre_path),
            "post_h5ad": str(post_path),
            "detections_parquet": str(parquet_path),
        },
    }
    runtime_path = task_root / "runtime" / f"harmony_real_{name}.json"
    runtime_path.parent.mkdir(parents=True, exist_ok=True)
    runtime_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(json.dumps(manifest, indent=2))
    logger.info("Done. Detections: %s", parquet_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
