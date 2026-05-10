"""
Unified integration-runs harness.

Run 9 integration methods on a list of AnnData datasets, saving:
    shared/integrations/<method>/<dataset>.h5ad   # X_integrated + obs/var
    shared/runtime/<method>_<dataset>.json        # wallclock + peak VRAM + RAM

Design:
    - One AnnData in → one AnnData out.
    - Every method is wrapped in `run(adata, batch_key, seed) -> (emb, elapsed_s, peak_vram_gb)`.
    - Methods self-select CPU or GPU. Concurrency is orchestrated at the caller level
      by pinning CUDA_VISIBLE_DEVICES (see experimental_plan_v1.md §4.2).
    - Deterministic under cfg.seed.

Methods implemented: Harmony, Seurat-RPCA (via scanpy + harmonypy fallback, or rpy2 bridge),
                     Scanorama, BBKNN, scVI, scANVI, scDML, scDREAMER, scCRAFT.

Stub placeholders live in `methods/<name>.py` so each can be iterated independently.
This top-level file only coordinates, so a missing method should not break the others.
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib
import json
import logging
import pathlib
import time
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)

DEFAULT_METHODS = [
    "harmony",
    "scanorama",
    "bbknn",
    "scvi",
    "scanvi",
    "scdml",
    "scdreamer",
    "sccraft",
    "seurat_rpca",
]


@dataclasses.dataclass
class HarnessConfig:
    batch_key: str = "batch"
    label_key: str | None = "cell_type"
    seed: int = 0
    out_root: str = "/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared"
    n_top_genes: int = 2000
    n_pcs: int = 50


def _peak_vram_gb() -> float:
    try:
        import torch
        if torch.cuda.is_available():
            return torch.cuda.max_memory_allocated() / (1024 ** 3)
    except ImportError:
        pass
    return 0.0


def _run_method(method: str, adata, cfg: HarnessConfig) -> tuple[np.ndarray, float, float]:
    mod = importlib.import_module(f"workspace.code.harness.methods.{method}")
    start = time.perf_counter()
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.reset_peak_memory_stats()
    except ImportError:
        pass
    emb = mod.run(adata, batch_key=cfg.batch_key, seed=cfg.seed)
    elapsed = time.perf_counter() - start
    return emb, elapsed, _peak_vram_gb()


def run_dataset(
    adata,
    dataset_name: str,
    methods: list[str],
    cfg: HarnessConfig,
) -> dict[str, Any]:
    """Run all methods on one dataset. Returns a dict of per-method outcomes."""
    results: dict[str, Any] = {}
    out_root = pathlib.Path(cfg.out_root)
    (out_root / "integrations").mkdir(parents=True, exist_ok=True)
    (out_root / "runtime").mkdir(parents=True, exist_ok=True)

    for method in methods:
        logger.info("Running %s on %s", method, dataset_name)
        try:
            emb, elapsed, vram = _run_method(method, adata, cfg)
            out = adata.copy()
            out.obsm["X_integrated"] = emb
            out_path = out_root / "integrations" / method / f"{dataset_name}.h5ad"
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out.write_h5ad(out_path)

            rt_path = out_root / "runtime" / f"{method}_{dataset_name}.json"
            rt_path.write_text(
                json.dumps(
                    {
                        "method": method,
                        "dataset": dataset_name,
                        "wallclock_s": elapsed,
                        "peak_vram_gb": vram,
                        "n_cells": int(adata.n_obs),
                        "n_genes": int(adata.n_vars),
                        "seed": cfg.seed,
                    },
                    indent=2,
                )
            )
            results[method] = {"emb_shape": emb.shape, "wallclock_s": elapsed, "vram_gb": vram}
        except Exception as exc:
            logger.exception("method %s on %s failed", method, dataset_name)
            results[method] = {"error": repr(exc)}

    return results


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--dataset", required=True)
    p.add_argument("--methods", nargs="+", default=DEFAULT_METHODS)
    p.add_argument("--batch-key", default="batch")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")
    import anndata as ad
    adata = ad.read_h5ad(args.h5ad)

    cfg = HarnessConfig(batch_key=args.batch_key, seed=args.seed)
    results = run_dataset(adata, args.dataset, args.methods, cfg)
    print(json.dumps({"dataset": args.dataset, "results": results}, indent=2))


if __name__ == "__main__":
    main()
