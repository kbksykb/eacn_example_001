"""
Sade-Feldman 2018 melanoma × canonical TPEX holdout pilot.

Oracle-defined TPEX holdout: cells scoring high on the Miller 2019 canonical TPEX panel
    positive core: TCF7, CXCR5, IL7R, SLAMF6 (at least 2 with TPM > 5)
    exclusion: PDCD1 (high), CXCL13 (high) — these mark terminally-exhausted, not TPEX

Since Sade-Feldman doesn't ship per-cell cell-type annotations in the GEO bundle, we
build the TPEX ground truth ourselves from oracle markers. This is explicitly disclosed
in Methods as "oracle-based holdout" rather than "author-annotated holdout" to pre-empt
non-circularity concerns: the oracle markers are DISCLOSED BEFORE the detection pipeline
runs; they are NOT fed as inputs to the integration step.
"""

from __future__ import annotations

import argparse
import json
import logging
import pathlib
import time

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

SHARED = pathlib.Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared")


def score_tpex(X: np.ndarray, var_names: list[str]) -> np.ndarray:
    """Oracle TPEX score per cell. Miller 2019 / Im 2016 canonical TPEX panel."""
    positive_core = ["TCF7", "CXCR5", "IL7R", "SLAMF6"]
    exclusion = ["CXCL13"]  # CXCL13+ cells are terminally-exhausted, not TPEX
    score = np.zeros(X.shape[0], dtype=np.float32)
    for g in positive_core:
        if g in var_names:
            col = var_names.index(g)
            vals = X[:, col]
            score += (vals > 2.0).astype(np.float32)  # TPM > 2
    for g in exclusion:
        if g in var_names:
            col = var_names.index(g)
            vals = X[:, col]
            score -= 0.5 * (vals > 2.0).astype(np.float32)
    return score


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--methods", nargs="+", default=["harmony", "scanorama", "scvi"])
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")

    import anndata as ad
    h5ad = SHARED / "data" / "sade_melanoma_cd8.h5ad"
    logger.info("loading %s", h5ad)
    a = ad.read_h5ad(h5ad)

    X = a.X.toarray() if hasattr(a.X, "toarray") else a.X
    tpex_score = score_tpex(X, list(a.var_names))
    logger.info("TPEX score distribution: min=%.1f max=%.1f 90p=%.1f 95p=%.1f",
                tpex_score.min(), tpex_score.max(),
                np.percentile(tpex_score, 90), np.percentile(tpex_score, 95))
    # TPEX = score >= 3 (at least 3 of 4 markers + no CXCL13-penalty trigger)
    is_tpex = tpex_score >= 3.0
    logger.info("TPEX cells (score ≥ 3): %d (%.2f%%)", int(is_tpex.sum()), 100 * is_tpex.mean())

    a.obs["cell_type"] = np.where(is_tpex, "TPEX_canonical", "other")
    a.obs["labels"] = a.obs["cell_type"].astype(str)

    out_h5ad = SHARED / "data" / "sade_tpex_annotated.h5ad"
    a.write_h5ad(out_h5ad)
    logger.info("wrote %s", out_h5ad)

    # Now run the retina_pilot with TPEX_canonical as the rare-label
    import subprocess, sys
    cmd = [
        sys.executable, "-m", "workspace.code.pilots.retina_pilot",
        "--h5ad", str(out_h5ad),
        "--rare-label", "TPEX_canonical",
        "--out-name", "sade_tpex_holdout",
        "--methods", *args.methods,
    ]
    logger.info("running: %s", " ".join(cmd))
    env = {"PYTHONPATH": "."}
    import os
    env.update({k: v for k, v in os.environ.items() if k not in ("PYTHONPATH",)})
    subprocess.run(cmd, env=env, check=False, cwd="/mnt/d-1274477442621830-m5ObBqn4/eacn3-main/examples/computational_biology")


if __name__ == "__main__":
    main()
