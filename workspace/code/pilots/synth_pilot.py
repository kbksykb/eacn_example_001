"""
End-to-end pilot on synthetic pancreas-like data.

Validates the complete pipeline plumbing before running on the real scIB pancreas:
    simulate → HVG + PCA → Leiden over-cluster (motif candidates) →
    Harmony + Scanorama + scVI → REAL channels (mutual-kNN purity +
    Procrustes displacement + bootstrap-stable) → per-motif L + LossRate@k →
    parquet output DATA_SPEC-compliant.

Runs fully in ≤ a few minutes on CPU/GPU. All outputs written to:
    /mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/data/pancreas_synth.h5ad      (pre)
    /mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/integrations/<method>/pancreas_synth.h5ad  (post)
    /mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/integrations/manifest.csv
    /mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/detections/real/<method>_pancreas_synth.parquet
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import pathlib
import resource
import subprocess
import time

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

SHARED = pathlib.Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared")


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(pathlib.Path(__file__).resolve().parents[4]), "rev-parse", "HEAD"],
            text=True,
        ).strip()[:12]
    except Exception:
        return "nogit"


def _peak_ram_gb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)


def _peak_vram_gb() -> float:
    try:
        import torch
        if torch.cuda.is_available():
            return torch.cuda.max_memory_allocated() / (1024 ** 3)
    except ImportError:
        pass
    return 0.0


def _reset_vram():
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.reset_peak_memory_stats()
    except ImportError:
        pass


def _pca_project(embedding: np.ndarray, d: int, seed: int) -> np.ndarray:
    """Project/trim an embedding to exactly d dimensions via PCA."""
    if embedding.shape[1] == d:
        return embedding
    if embedding.shape[1] < d:
        pad = np.zeros((embedding.shape[0], d - embedding.shape[1]), dtype=embedding.dtype)
        return np.concatenate([embedding, pad], axis=1)
    from sklearn.decomposition import PCA
    pca = PCA(n_components=d, random_state=seed)
    return pca.fit_transform(embedding).astype(np.float32)


def preprocess(adata, batch_key: str, seed: int):
    import scanpy as sc

    a = adata.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=1000, batch_key=batch_key, flavor="seurat")
    a.var["highly_variable"] = a.var["highly_variable"].fillna(False)
    hvg = a[:, a.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=50, random_state=seed)
    a.obsm["X_uncorrected_pca"] = _pca_project(hvg.obsm["X_pca"], 50, seed)
    # Downstream sc.pp.neighbors needs X_pca alias
    a.obsm["X_pca"] = a.obsm["X_uncorrected_pca"]
    return a


def over_cluster(adata, seed: int, resolution: float = 3.0) -> np.ndarray:
    import scanpy as sc

    sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca", random_state=seed)
    sc.tl.leiden(
        adata,
        resolution=resolution,
        random_state=seed,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    return adata.obs["leiden"].astype(int).to_numpy()


def run_harmony(adata, batch_key: str, seed: int) -> tuple[np.ndarray, dict]:
    import harmonypy as hm

    ho = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        vars_use=[batch_key],
        random_state=seed,
    )
    # harmonypy's Z_corr is (d, n); transpose to (n, d).
    z = np.asarray(ho.Z_corr)
    n = adata.n_obs
    if z.shape[0] == n:
        emb = z.astype(np.float32)
    elif z.shape[1] == n:
        emb = z.T.astype(np.float32)
    else:
        raise RuntimeError(f"harmonypy Z_corr shape {z.shape} does not match n_obs {n}")
    return _pca_project(emb, 32, seed), {
        "method": "harmony",
        "method_version": getattr(hm, "__version__", "unknown"),
    }


def run_scanorama(adata, batch_key: str, seed: int) -> tuple[np.ndarray, dict]:
    import scanorama

    batches = adata.obs[batch_key].astype(str)
    uniq = sorted(batches.unique())
    splits = [adata[batches == b].copy() for b in uniq]
    out = scanorama.correct_scanpy(splits, return_dimred=True)
    # Different Scanorama versions return either (corrected, dimred) or just corrected.
    if isinstance(out, tuple) and len(out) == 2:
        corrected = out[0]
    else:
        corrected = out
    emb = np.empty((adata.n_obs, corrected[0].obsm["X_scanorama"].shape[1]), dtype=np.float32)
    for b, piece in zip(uniq, corrected):
        emb[np.where(batches == b)[0]] = piece.obsm["X_scanorama"]
    return _pca_project(emb, 32, seed), {
        "method": "scanorama",
        "method_version": getattr(scanorama, "__version__", "unknown"),
    }


def run_scvi(adata, batch_key: str, seed: int) -> tuple[np.ndarray, dict]:
    import scvi

    scvi.settings.seed = seed
    a = adata.copy()
    import scanpy as sc
    sc.pp.highly_variable_genes(a, n_top_genes=1000, batch_key=batch_key, flavor="seurat_v3")
    a = a[:, a.var["highly_variable"]].copy()
    scvi.model.SCVI.setup_anndata(a, batch_key=batch_key)
    model = scvi.model.SCVI(a, n_latent=32)
    model.train(max_epochs=80, early_stopping=True, check_val_every_n_epoch=10, enable_progress_bar=False)
    emb = model.get_latent_representation().astype(np.float32)
    return emb, {"method": "scvi", "method_version": scvi.__version__}


METHODS = {
    "harmony": run_harmony,
    "scanorama": run_scanorama,
    "scvi": run_scvi,
}


def write_pre(adata, dataset_name: str, run_info: dict) -> pathlib.Path:
    SHARED.mkdir(parents=True, exist_ok=True)
    data_dir = SHARED / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    out = data_dir / f"{dataset_name}.h5ad"
    a = adata.copy()
    # DATA_SPEC required fields if not yet present
    if "cell_type_public" not in a.obs.columns:
        a.obs["cell_type_public"] = "NA"
    a.uns["reproducibility"] = run_info
    a.write_h5ad(out)
    logger.info("pre-integration written → %s", out)
    return out


def write_post(adata_pre, emb_int: np.ndarray, method: str, dataset_name: str, run_info: dict):
    a = adata_pre.copy()
    a.obsm["X_integrated"] = emb_int
    a.uns["method"] = run_info.get("method", method)
    a.uns["method_version"] = run_info.get("method_version", "unknown")
    a.uns["method_config"] = run_info.get("method_config", {})
    a.uns["runtime"] = run_info["runtime"]
    out_dir = SHARED / "integrations" / method
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{dataset_name}.h5ad"
    a.write_h5ad(out)
    logger.info("post-integration written → %s", out)
    return out


def run_real_channels(adata_pre, emb_int: np.ndarray, motif_candidates: np.ndarray):
    import pandas as pd
    from workspace.code.real_channels import ChannelConfig, loss_rate_at_k, score_purity_and_procrustes

    batch_int = pd.Categorical(adata_pre.obs["batch"]).codes.astype(np.int32)
    d = emb_int.shape[1]
    pre_full = adata_pre.obsm["X_uncorrected_pca"].astype(np.float32)
    # Truncate pre to same dim as post so Procrustes is well-posed.
    pre_trunc = pre_full[:, :d] if pre_full.shape[1] >= d else _pca_project(pre_full, d, 0)
    cfg = ChannelConfig(k=15, n_bootstrap=30, rare_abundance_threshold=0.10, seed=0)
    reports = score_purity_and_procrustes(
        emb_pre=pre_trunc,
        emb_post=emb_int.astype(np.float32),
        candidate_labels=motif_candidates,
        batch_labels=batch_int,
        cfg=cfg,
    )
    return reports, loss_rate_at_k(reports, k=1), loss_rate_at_k(reports, k=3)


def ground_truth_rare_candidates(motif_candidates: np.ndarray, hidden_labels: np.ndarray, rare_label: str, min_overlap: int = 3) -> set[int]:
    is_rare = np.asarray([str(x).lower().strip() == rare_label.lower().strip() for x in hidden_labels])
    out = set()
    for cid in np.unique(motif_candidates):
        mask = motif_candidates == cid
        if mask.sum() == 0:
            continue
        if is_rare[mask].sum() >= min_overlap and is_rare[mask].sum() / mask.sum() > 0.5:
            out.add(int(cid))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--methods", nargs="+", default=["harmony", "scanorama", "scvi"])
    p.add_argument("--dataset", default="pancreas_synth")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--rare-label", default="epsilon")
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")

    from workspace.code.datasets.synthetic import SimConfig, simulate

    sim_cfg = SimConfig(n_cells_per_batch=2000, n_batches=3, seed=args.seed)
    adata_raw = simulate(sim_cfg)
    hidden = adata_raw.obs["__hidden_label__"].astype(str).to_numpy()
    logger.info("rare cells in truth: %d", int(adata_raw.obs["heldout_rare"].sum()))

    logger.info("preprocessing…")
    adata_pre = preprocess(adata_raw, batch_key="batch", seed=args.seed)
    logger.info("over-clustering…")
    candidates = over_cluster(adata_pre, seed=args.seed)
    adata_pre.obs["motif_id_pre"] = candidates

    run_info = {
        "random_state": args.seed,
        "commit_sha": _git_sha(),
        "pip_freeze_hash": hashlib.sha256(b"scenv").hexdigest()[:16],  # placeholder
        "run_date": time.strftime("%Y-%m-%d"),
    }
    write_pre(adata_pre, args.dataset, run_info)

    manifest_rows = []
    summary_rows = []
    det_dir = SHARED / "detections" / "real"
    det_dir.mkdir(parents=True, exist_ok=True)

    for method in args.methods:
        if method not in METHODS:
            logger.warning("method %s not registered — skipping", method)
            continue
        logger.info("=== %s ===", method)
        _reset_vram()
        t0 = time.perf_counter()
        try:
            emb_int, method_meta = METHODS[method](adata_pre, batch_key="batch", seed=args.seed)
        except Exception:
            logger.exception("method %s failed", method)
            continue
        elapsed = time.perf_counter() - t0
        vram = _peak_vram_gb()
        ram = _peak_ram_gb()

        post_run_info = {
            "method": method_meta["method"],
            "method_version": method_meta["method_version"],
            "method_config": {"seed": args.seed, "d_projection": 32},
            "runtime": {
                "wall_clock_seconds": round(elapsed, 3),
                "peak_vram_gb": round(vram, 3),
                "peak_ram_gb": round(ram, 3),
                "n_cells_used": int(adata_pre.n_obs),
            },
        }
        write_post(adata_pre, emb_int, method, args.dataset, post_run_info)

        reports, lr1, lr3 = run_real_channels(adata_pre, emb_int, candidates)

        rare_ids = ground_truth_rare_candidates(candidates, hidden, args.rare_label)
        rows = []
        for r in reports:
            rows.append({
                "method": method,
                "dataset": args.dataset,
                "candidate_id": r.candidate_id,
                "abundance": r.abundance,
                "channel_mknn_purity_pre": r.purity_pre,
                "channel_mknn_purity_post": r.purity_post,
                "channel_mknn_nn_pre": r.nn_identity_pre,
                "channel_mknn_nn_post": r.nn_identity_post,
                "channel_proc_displacement": r.procrustes_displacement,
                "channel_boot_stable": r.bootstrap_fraction_stable,
                "loss_probability": r.loss_probability,
                "bootstrap_ci_low": float(np.nan),
                "bootstrap_ci_high": float(np.nan),
                "ground_truth_label": (r.candidate_id in rare_ids),
                "seed": args.seed,
                "n_cells": int(adata_pre.n_obs),
            })
        df = pd.DataFrame(rows)
        parquet_path = det_dir / f"{method}_{args.dataset}.parquet"
        df.to_parquet(parquet_path)
        logger.info(
            "%s: wrote %s (%d rows, %d rare-truth rows), LossRate@1=%.3f, @3=%.3f, wall=%.1fs, vram=%.1fGB",
            method, parquet_path, len(df), int(df["ground_truth_label"].sum()), lr1, lr3, elapsed, vram,
        )

        manifest_rows.append({
            "method": method,
            "dataset": args.dataset,
            "wall_clock_seconds": round(elapsed, 3),
            "peak_vram_gb": round(vram, 3),
            "peak_ram_gb": round(ram, 3),
            "n_cells_used": int(adata_pre.n_obs),
            "package_version": method_meta["method_version"],
            "commit_sha": run_info["commit_sha"],
            "run_date": run_info["run_date"],
        })
        summary_rows.append({
            "method": method,
            "LossRate@1": lr1,
            "LossRate@3": lr3,
            "n_candidates": len(reports),
            "n_rare_truth_candidates": len(rare_ids),
            "max_loss_prob_on_rare_truth": float(df.loc[df["ground_truth_label"], "loss_probability"].max() if df["ground_truth_label"].any() else 0.0),
            "max_loss_prob_on_non_rare_truth": float(df.loc[~df["ground_truth_label"], "loss_probability"].max() if (~df["ground_truth_label"]).any() else 0.0),
        })

    # Manifest
    manifest_path = SHARED / "integrations" / "manifest.csv"
    mdf = pd.DataFrame(manifest_rows)
    if manifest_path.exists():
        existing = pd.read_csv(manifest_path)
        existing = existing[~((existing["method"].isin(mdf["method"])) & (existing["dataset"] == args.dataset))]
        mdf = pd.concat([existing, mdf], ignore_index=True)
    mdf.to_csv(manifest_path, index=False)
    logger.info("manifest → %s", manifest_path)

    # Summary
    print("\n=== SUMMARY ===")
    sdf = pd.DataFrame(summary_rows)
    print(sdf.to_string(index=False))
    results_dir = pathlib.Path(__file__).resolve().parents[3] / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    sdf.to_csv(results_dir / f"{args.dataset}_pilot_summary.csv", index=False)


if __name__ == "__main__":
    main()
