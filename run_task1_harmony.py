#!/usr/bin/env python3
"""
Task 1 — Harmony batch integration (harness benchmark spec).

Wraps workspace.code.harness with 2000 HVG, PCA 50, harmonypy on X_pca.
Outputs live under <out-root>/task1/ (separate from task 2).

Usage (from repo root):
    python run_task1_harmony.py --input workspace/Galaxy1.h5ad
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time

import anndata as ad

from pipeline_common import (
    add_common_input_args,
    dataset_name_from_input,
    git_sha_short,
    resolve_input_path,
    resolve_out_root,
    validate_batch_key,
)
from workspace.code.harness import HarnessConfig, run_dataset

logger = logging.getLogger(__name__)


def main() -> int:
    p = argparse.ArgumentParser(
        description="Task 1: HVG + PCA + Harmony via harness (benchmark integration spec)."
    )
    add_common_input_args(p)
    args = p.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    input_path = resolve_input_path(args.input)
    name = dataset_name_from_input(input_path, args.name)
    out_root = resolve_out_root(args.out_root)
    task_root = out_root / "task1"

    logger.info("Task 1 — Harmony (harness)")
    logger.info("  input: %s", input_path)
    logger.info("  name:  %s", name)
    logger.info("  out:   %s", task_root)

    adata = ad.read_h5ad(input_path)
    validate_batch_key(adata, args.batch_key)

    t0 = time.perf_counter()
    cfg = HarnessConfig(
        batch_key=args.batch_key,
        seed=args.seed,
        out_root=str(task_root),
    )
    results = run_dataset(adata, name, ["harmony"], cfg)
    elapsed = time.perf_counter() - t0

    if "error" in results.get("harmony", {}):
        logger.error("Harmony failed: %s", results["harmony"]["error"])
        return 1

    h5ad_out = task_root / "integrations" / "harmony" / f"{name}.h5ad"
    runtime_out = task_root / "runtime" / f"harmony_{name}.json"

    manifest = {
        "task": "task1_harmony",
        "pipeline": "harness/methods/harmony.py",
        "spec": "2000 HVG, PCA 50, X_integrated same dimension as PCA",
        "input_h5ad": str(input_path),
        "dataset_name": name,
        "batch_key": args.batch_key,
        "seed": args.seed,
        "git_commit": git_sha_short(),
        "total_wallclock_s": elapsed,
        "harmony_result": results.get("harmony"),
        "outputs": {
            "integrated_h5ad": str(h5ad_out),
            "runtime_json": str(runtime_out),
        },
    }
    manifest_path = task_root / "runtime" / f"task1_{name}_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(json.dumps(manifest, indent=2))
    logger.info("Done. Integrated h5ad: %s", h5ad_out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
