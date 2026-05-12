"""
kang2024_melanoma_slice.py — Melanoma slice integration + TPEX detection for Fig 5.

Merges 4 Kang 2024 melanoma cohorts:
  - mel_GSE200218 (Sade-Feldman 2018 — primary melanoma, ICI response)
  - mel_uv_GSE138433 (uveal melanoma)
  - mel_uv_GSE139829 (uveal melanoma)
  - mel_uv_GSE176029 (uveal melanoma)

Protocol:
  1. Decompress + load each h5ad.
  2. Merge, harmonize obs columns.
  3. Preprocess (HVG, PCA).
  4. Run Harmony integration.
  5. Detect TPEX (TCF7/TCF1+ CD8 T cells) via OT channel.
  6. Detect LAMP3+ mregDC via OT channel.
  7. Save integrated h5ad + detection parquet for Immunology annotation.

Output:
  workspace/results/kang2024_melanoma/melanoma_harmony.h5ad
  workspace/results/kang2024_melanoma/melanoma_detection.parquet
  workspace/results/kang2024_melanoma/melanoma_summary.json
"""
from __future__ import annotations

import json
import lzma
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

ATLAS_DIR = Path("/mnt/kang_2024_pan_cancer_atlas/extracted/atlas_dataset")
WORKSPACE = Path(__file__).parent.parent.parent
RESULTS_DIR = WORKSPACE / "results" / "kang2024_melanoma"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

MEL_FILES = [
    "mel_GSE200218.h5ad.xz",
    "mel_uv_GSE138433.h5ad.xz",
    "mel_uv_GSE139829.h5ad.xz",
    "mel_uv_GSE176029.h5ad.xz",
]

# TPEX markers (Miller 2019, canonical TCF1+ exhausted CD8)
TPEX_MARKERS = ["TCF7", "TOX", "PDCD1", "HAVCR2", "ENTPD1", "CXCR5"]
# LAMP3+ mregDC markers
LAMP3_MARKERS = ["LAMP3", "CCR7", "FSCN1", "MARCKSL1", "BIRC3"]

BATCH_KEY_CANDIDATES = ["batch", "sample", "Sample", "donor_id", "Patient",
                        "patient", "Patient_Organ_Tissue"]


def _find_batch_key(adata) -> str:
    for k in BATCH_KEY_CANDIDATES:
        if k in adata.obs.columns:
            return k
    for k in adata.obs.columns:
        if adata.obs[k].dtype.name in ("category", "object"):
            return k
    raise ValueError(f"No batch key. obs: {list(adata.obs.columns)}")


def _decompress_load(xz_path: Path) -> sc.AnnData:
    """Decompress .xz on-the-fly and load h5ad."""
    h5_path = RESULTS_DIR / xz_path.stem  # strip .xz
    if not h5_path.exists():
        print(f"  Decompressing {xz_path.name}...")
        with lzma.open(xz_path) as f_in, open(h5_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    return sc.read_h5ad(h5_path)


def _json_safe(obj):
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
    try:
        import torch
        if isinstance(obj, (torch.dtype, torch.device)):
            return str(obj)
    except ImportError:
        pass
    try:
        json.dumps(obj)
        return obj
    except (TypeError, ValueError):
        return str(obj)


def load_and_merge() -> sc.AnnData:
    """Load all melanoma cohorts and merge."""
    adatas = []
    for fname in MEL_FILES:
        fpath = ATLAS_DIR / fname
        if not fpath.exists():
            print(f"  WARNING: {fname} not found, skipping")
            continue
        print(f"  Loading {fname}...")
        adata = _decompress_load(fpath)
        adata.obs["source_file"] = fname.replace(".h5ad.xz", "")
        # Standardize batch key
        bk = _find_batch_key(adata)
        adata.obs["batch"] = adata.obs[bk].astype(str) + "_" + fname[:10]
        print(f"    n_obs={adata.n_obs}, batch_key={bk}, n_batches={adata.obs[bk].nunique()}")
        adatas.append(adata)

    if not adatas:
        raise RuntimeError("No melanoma files loaded")

    print(f"  Merging {len(adatas)} cohorts...")
    merged = sc.concat(adatas, join="inner", label="cohort",
                       keys=[a.obs["source_file"].iloc[0] for a in adatas])
    print(f"  Merged: n_obs={merged.n_obs}, n_vars={merged.n_vars}")
    return merged


def build_candidate_labels(adata, markers: list[str], top_pct: float = 0.01) -> np.ndarray:
    """Score cells by marker co-expression; top top_pct = candidate."""
    score = np.zeros(adata.n_obs)
    found = []
    for gene in markers:
        if gene in adata.var_names:
            idx = list(adata.var_names).index(gene)
            score += np.asarray(adata.X[:, idx]).ravel()
            found.append(gene)
    if not found:
        print(f"    WARNING: none of {markers} found in var_names")
        return np.full(adata.n_obs, -1, dtype=int)
    threshold = np.percentile(score, 100 * (1 - top_pct))
    labels = np.where(score >= threshold, 0, -1).astype(int)
    n_pos = int((labels == 0).sum())
    print(f"    Markers found: {found}; n_candidate={n_pos} ({100*n_pos/adata.n_obs:.2f}%)")
    return labels


def run_ot_detection(adata, z_pre, z_post, candidate_labels, label_name: str) -> dict:
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from real_channels.ot_channel import OTChannelConfig, score_ot

    n_pos = int((candidate_labels == 0).sum())
    if n_pos < 10:
        return {"skipped": True, "reason": f"too few {label_name} cells: {n_pos}"}

    batch_labels = adata.obs["batch"].astype("category").cat.codes.values

    # Use CPU to avoid O(n²) GPU OOM for n>10k; sklearn kNN is used internally
    # when device="cpu" and n > local_density threshold.
    # Subsample to 15k for tractable CPU pairwise distance if needed.
    n = z_pre.shape[0]
    if n > 15_000:
        rng = np.random.default_rng(42)
        idx = rng.choice(n, size=15_000, replace=False)
        # Keep all candidate cells + random sample of non-candidates
        cand_idx = np.where(candidate_labels >= 0)[0]
        non_cand_idx = np.where(candidate_labels < 0)[0]
        n_fill = max(0, 15_000 - len(cand_idx))
        if n_fill > 0 and len(non_cand_idx) > 0:
            fill_idx = rng.choice(non_cand_idx,
                                  size=min(n_fill, len(non_cand_idx)),
                                  replace=False)
            idx = np.concatenate([cand_idx, fill_idx])
        else:
            idx = cand_idx
        z_pre_sub = z_pre[idx]
        z_post_sub = z_post[idx]
        cand_sub = candidate_labels[idx]
        batch_sub = batch_labels[idx]
        print(f"    Subsampled to {len(idx)} cells for OT detection")
    else:
        z_pre_sub, z_post_sub = z_pre, z_post
        cand_sub, batch_sub = candidate_labels, batch_labels

    cfg = OTChannelConfig(
        n_permutations=200,
        k_neighbors=30,
        device="cpu",  # avoid GPU OOM for large n
        seed=42,
    )
    result = score_ot(z_pre_sub, z_post_sub, cand_sub, batch_sub, cfg)
    result["label_name"] = label_name
    result["n_candidate"] = n_pos
    result["n_total"] = adata.n_obs
    return result


def _gpu_available() -> bool:
    try:
        import torch
        return torch.cuda.is_available()
    except ImportError:
        return False


def main():
    t_start = time.time()
    print("=== Kang 2024 Melanoma Slice ===")

    adata = load_and_merge()

    # Preprocess
    print("  Preprocessing...")
    if adata.X.max() > 100:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
    z_pre = adata.obsm["X_pca"].copy()
    print(f"  After HVG: n_vars={adata.n_vars}")

    # Harmony
    import harmonypy as hm
    print(f"  Running Harmony (n_batches={adata.obs['batch'].nunique()})...")
    t_h = time.time()
    ho = hm.run_harmony(z_pre, adata.obs, "batch", max_iter_harmony=20,
                        random_state=0)
    z_post = np.array(ho.Z_corr)
    print(f"  Harmony done in {time.time()-t_h:.1f}s")

    adata.obsm["X_harmony"] = z_post

    # Detection
    print("  Detecting TPEX...")
    tpex_labels = build_candidate_labels(adata, TPEX_MARKERS, top_pct=0.02)
    tpex_result = run_ot_detection(adata, z_pre, z_post, tpex_labels, "TPEX")

    print("  Detecting LAMP3+ mregDC...")
    lamp3_labels = build_candidate_labels(adata, LAMP3_MARKERS, top_pct=0.01)
    lamp3_result = run_ot_detection(adata, z_pre, z_post, lamp3_labels, "LAMP3_mregDC")

    t_total = time.time() - t_start

    # Save detection parquet
    rows = []
    for res, name in [(tpex_result, "TPEX"), (lamp3_result, "LAMP3_mregDC")]:
        if res.get("skipped"):
            continue
        for i, cid in enumerate(res["candidate_ids"]):
            rows.append({
                "motif_id": f"{name}_{cid}",
                "label_name": name,
                "n_cells": res["n_candidate"],
                "n_total": res["n_total"],
                "ot_stat": float(res["stat"][i]),
                "ot_p_value": float(res["p_values"][i]),
                "ci_low": float(res["ci_low"][i]),
                "ci_high": float(res["ci_high"][i]),
                "loss_probability": 1.0 - float(res["p_values"][i]),
            })
    if rows:
        det_df = pd.DataFrame(rows)
        det_df.to_parquet(RESULTS_DIR / "melanoma_detection.parquet", index=False)
        print(f"  Detection parquet: {len(rows)} rows")

    # Save summary
    summary = {
        "n_obs": adata.n_obs,
        "n_vars_hvg": adata.n_vars,
        "n_batches": adata.obs["batch"].nunique(),
        "cohorts": MEL_FILES,
        "t_total_s": round(t_total, 1),
        "tpex_detection": _json_safe(tpex_result),
        "lamp3_detection": _json_safe(lamp3_result),
    }
    with open(RESULTS_DIR / "melanoma_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Save integrated h5ad (without raw counts to save space)
    print("  Saving integrated h5ad...")
    adata.write_h5ad(RESULTS_DIR / "melanoma_harmony.h5ad", compression="gzip")

    print(f"  DONE: {adata.n_obs} cells, {t_total:.1f}s total")
    print(f"  TPEX p={tpex_result.get('p_values', [None])[0]}")
    print(f"  LAMP3 p={lamp3_result.get('p_values', [None])[0]}")


if __name__ == "__main__":
    main()
