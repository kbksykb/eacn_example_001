"""Shared helpers for run_task1_harmony.py and run_task2_harmony_real.py."""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent

# Recommended place to put input .h5ad files (any path works via --input).
DEFAULT_INPUT_DIR = REPO_ROOT / "workspace" / "data" / "input"

# Server shared root from experimental_plan_v1.md (used only as README / env hint).
SERVER_SHARED = Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared")


def resolve_out_root(explicit: str | None) -> Path:
    if explicit:
        return Path(explicit).expanduser().resolve()
    env = os.environ.get("EACN_OUT_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    return (REPO_ROOT / "output").resolve()


def resolve_input_path(path: str) -> Path:
    p = Path(path).expanduser()
    if not p.is_absolute():
        p = (REPO_ROOT / p).resolve()
    if not p.exists():
        raise FileNotFoundError(f"Input h5ad not found: {p}")
    if p.suffix.lower() not in {".h5ad", ".h5"}:
        raise ValueError(f"Expected an .h5ad file, got: {p}")
    return p


def dataset_name_from_input(input_path: Path, explicit: str | None) -> str:
    if explicit:
        return explicit
    return input_path.stem


def validate_batch_key(adata, batch_key: str) -> None:
    if batch_key not in adata.obs.columns:
        cols = ", ".join(map(str, adata.obs.columns[:20]))
        suffix = " …" if len(adata.obs.columns) > 20 else ""
        raise KeyError(
            f"Batch column {batch_key!r} not in adata.obs. "
            f"Available columns: {cols}{suffix}. Use --batch-key."
        )


def git_sha_short() -> str:
    try:
        return (
            subprocess.check_output(
                ["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"],
                text=True,
            )
            .strip()[:12]
        )
    except Exception:
        return "nogit"


def add_common_input_args(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "--input",
        required=True,
        help="Path to input .h5ad (absolute, or relative to repo root). "
        f"Suggested folder: {DEFAULT_INPUT_DIR.relative_to(REPO_ROOT)}/",
    )
    p.add_argument(
        "--name",
        default=None,
        help="Dataset name for output filenames (default: input file stem, e.g. Galaxy1).",
    )
    p.add_argument(
        "--batch-key",
        default="batch",
        help="obs column with batch labels (default: batch).",
    )
    p.add_argument(
        "--out-root",
        default=None,
        help="Root output directory (default: $EACN_OUT_ROOT or ./output/). "
        "Task outputs go under <out-root>/task1/ and <out-root>/task2/.",
    )
    p.add_argument("--seed", type=int, default=0, help="Random seed (default: 0).")
