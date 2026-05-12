"""
kang2024_fig4_pipeline.py — Kang 2024 pan-cancer atlas scalability run for Fig 4.

Protocol:
  1. Load each NMF compartment h5ad (T/NK, myeloid, B, mesenchymal, epithelial).
  2. Preprocess: HVG selection, PCA.
  3. Run Harmony integration (GPU-accelerated via harmonypy).
  4. Apply OT-channel detection on LAMP3+ mregDC holdout within myeloid compartment.
  5. Report: cell counts, runtime, peak VRAM, detection p-values.

GPU allocation: one compartment per GPU (0-4), all launched in parallel.
Output: workspace/results/kang2024_fig4_scalability.md + per-compartment parquets.

Usage:
    CUDA_VISIBLE_DEVICES=0 python kang2024_fig4_pipeline.py --compartment tnk
    CUDA_VISIBLE_DEVICES=1 python kang2024_fig4_pipeline.py --compartment myl
    # etc. — or use launch_all_compartments.sh to fire all 5 in parallel.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
NMF_DIR = Path("/mnt/kang_2024_pan_cancer_atlas/extracted/NMF_h5ad_decompressed")
WORKSPACE = Path(__file__).parent.parent.parent  # workspace/
RESULTS_DIR = WORKSPACE / "results" / "kang2024_fig4"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

COMPARTMENT_FILES = {
    "tnk":          "tnk_NMF.h5ad",
    "myl":          "myl_NMF.h5ad",
    "b":            "b_NMF.h5ad",
    "mesenchymal":  "mesenchymal_NMF.h5ad",
    "epi":          "epi_NMF.h5ad",
}

# LAMP3+ mregDC marker genes (from Cheng 2021 + Immunology agent catalog)
LAMP3_MARKERS = ["LAMP3", "CCR7", "FSCN1", "MARCKSL1", "BIRC3"]

# Batch key to use for integration (Kang atlas uses 'batch' or 'sample')
BATCH_KEY_CANDIDATES = ["batch", "sample", "Sample", "donor_id", "study"]


def _find_batch_key(adata) -> str:
    for k in BATCH_KEY_CANDIDATES:
        if k in adata.obs.columns:
            return k
    # Fall back to first categorical column
    for k in adata.obs.columns:
        if adata.obs[k].dtype.name in ("category", "object"):
            return k
    raise ValueError(f"No batch key found. obs columns: {list(adata.obs.columns)}")


def _peak_vram_mb() -> float:
    try:
        import torch
        if torch.cuda.is_available():
            return torch.cuda.max_memory_allocated() / 1e6
    except ImportError:
        pass
    return 0.0


def preprocess(adata, n_hvg: int = 3000, n_pcs: int = 50) -> sc.AnnData:
    """Standard preprocessing: normalize, log1p, HVG, PCA."""
    t0 = time.time()
    print(f"  [preprocess] n_obs={adata.n_obs}, n_vars={adata.n_vars}")

    # Only normalize if raw counts present
    if adata.X.max() > 100:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, batch_key=None, subset=True)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")
    print(f"  [preprocess] done in {time.time()-t0:.1f}s")
    return adata


def run_harmony(adata, batch_key: str, n_pcs: int = 50) -> np.ndarray:
    """Run Harmony on PCA embedding. Returns corrected embedding (n_cells, n_pcs)."""
    import harmonypy as hm
    t0 = time.time()
    print(f"  [harmony] batch_key={batch_key}, n_batches={adata.obs[batch_key].nunique()}")
    ho = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        batch_key,
        max_iter_harmony=20,
        random_state=0,
    )
    # harmonypy 2.x returns Z_corr as (n_cells, n_pcs)
    z_post = np.array(ho.Z_corr)  # (n_cells, n_pcs)
    print(f"  [harmony] done in {time.time()-t0:.1f}s")
    return z_post


def detect_lamp3_holdout(adata, z_pre: np.ndarray, z_post: np.ndarray) -> dict:
    """Run OT-channel detection on LAMP3+ mregDC holdout within myeloid compartment."""
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from real_channels.ot_channel import OTChannelConfig, score_ot

    # Build candidate labels: LAMP3+ cells = candidate 0, rest = -1
    lamp3_score = np.zeros(adata.n_obs)
    for gene in LAMP3_MARKERS:
        if gene in adata.var_names:
            idx = adata.var_names.get_loc(gene)
            lamp3_score += np.asarray(adata.X[:, idx]).ravel()

    # Top 1% by LAMP3 score = rare holdout
    threshold = np.percentile(lamp3_score, 99)
    candidate_labels = np.where(lamp3_score >= threshold, 0, -1)
    n_lamp3 = int((candidate_labels == 0).sum())
    print(f"  [lamp3_holdout] n_lamp3={n_lamp3} ({100*n_lamp3/adata.n_obs:.2f}%)")

    if n_lamp3 < 10:
        return {"skipped": True, "reason": "too few LAMP3+ cells"}

    batch_key = _find_batch_key(adata)
    batch_labels = adata.obs[batch_key].astype("category").cat.codes.values

    cfg = OTChannelConfig(
        n_permutations=200,
        k_neighbors=30,
        device="cuda" if _gpu_available() else "cpu",
        seed=42,
    )
    result = score_ot(z_pre, z_post, candidate_labels, batch_labels, cfg)
    result["n_lamp3"] = n_lamp3
    result["n_total"] = adata.n_obs
    result["frac_lamp3"] = n_lamp3 / adata.n_obs
    return result


def _gpu_available() -> bool:
    try:
        import torch
        return torch.cuda.is_available()
    except ImportError:
        return False


def _json_safe(obj):
    """Recursively convert numpy/torch types for JSON serialization."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    # Handle torch types (dtype, device, etc.)
    try:
        import torch
        if isinstance(obj, torch.dtype):
            return str(obj)
        if isinstance(obj, torch.device):
            return str(obj)
    except ImportError:
        pass
    # Fallback: convert to string for any other non-serializable type
    try:
        import json as _json
        _json.dumps(obj)
        return obj
    except (TypeError, ValueError):
        return str(obj)


def run_compartment(compartment: str, n_hvg: int = 3000, n_pcs: int = 50) -> dict:
    """Full pipeline for one compartment. Returns summary dict."""
    h5ad_path = NMF_DIR / COMPARTMENT_FILES[compartment]
    if not h5ad_path.exists():
        # Try .xz version (not yet decompressed)
        xz_path = Path(str(h5ad_path) + ".xz")
        if xz_path.exists():
            print(f"  [load] decompressing {xz_path.name} on-the-fly...")
            import lzma, shutil
            with lzma.open(xz_path) as f_in, open(h5ad_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        else:
            raise FileNotFoundError(f"Neither {h5ad_path} nor {xz_path} found")

    print(f"\n=== Compartment: {compartment} ===")
    t_start = time.time()

    print(f"  [load] reading {h5ad_path.name}...")
    adata = sc.read_h5ad(h5ad_path)
    print(f"  [load] n_obs={adata.n_obs}, n_vars={adata.n_vars}")
    print(f"  [load] obs columns: {list(adata.obs.columns)[:10]}")

    batch_key = _find_batch_key(adata)
    z_pre = adata.obsm.get("X_pca", None)

    adata = preprocess(adata, n_hvg=n_hvg, n_pcs=n_pcs)
    z_pre_pca = adata.obsm["X_pca"].copy()

    if z_pre is None:
        z_pre = z_pre_pca

    z_post = run_harmony(adata, batch_key=batch_key, n_pcs=n_pcs)
    adata.obsm["X_harmony"] = z_post

    t_integration = time.time() - t_start
    peak_vram = _peak_vram_mb()

    # Detection only on myeloid (has LAMP3+ mregDC)
    detection = {}
    if compartment == "myl":
        detection = detect_lamp3_holdout(adata, z_pre_pca, z_post)

    t_total = time.time() - t_start

    summary = {
        "compartment": compartment,
        "n_obs": adata.n_obs,
        "n_vars_hvg": adata.n_vars,
        "n_batches": adata.obs[batch_key].nunique(),
        "batch_key": batch_key,
        "t_integration_s": round(t_integration, 1),
        "t_total_s": round(t_total, 1),
        "peak_vram_mb": round(peak_vram, 1),
        "detection": detection,
    }

    # Save per-compartment result
    out_json = RESULTS_DIR / f"{compartment}_summary.json"
    with open(out_json, "w") as f:
        json.dump(_json_safe(summary), f, indent=2)

    # Save harmony embedding
    emb_out = RESULTS_DIR / f"{compartment}_harmony_embedding.npy"
    np.save(emb_out, z_post)

    print(f"  [done] {compartment}: {adata.n_obs} cells, {t_total:.1f}s total, "
          f"{peak_vram:.0f} MB VRAM")
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--compartment", choices=list(COMPARTMENT_FILES.keys()),
                        required=True)
    parser.add_argument("--n_hvg", type=int, default=3000)
    parser.add_argument("--n_pcs", type=int, default=50)
    args = parser.parse_args()

    summary = run_compartment(args.compartment, n_hvg=args.n_hvg, n_pcs=args.n_pcs)
    print(json.dumps(_json_safe(summary), indent=2))


if __name__ == "__main__":
    main()
