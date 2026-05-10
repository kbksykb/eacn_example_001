"""
Synthetic sweep across (abundance × overlap-strength) per Math's detectability envelope.

Fills the factorial grid BioSci asked for (Figure 3c-d in the outline):

    abundance π ∈ {0.005, 0.01, 0.02, 0.05} — 0.5% / 1% / 2% / 5%
    overlap   α ∈ {0.0, 0.2, 0.4, 0.6}       — epsilon markers overlap alpha markers at α strength
    methods   = {Harmony, Scanorama, scVI}

For each cell, emit per-method Fisher-fused p-values for the rare motif vs abundant motifs.

Outputs:
    workspace/results/synth_sweep/sweep_summary.csv (one row per abundance × α × method × seed)
    workspace/results/synth_sweep/sweep_detections.parquet (per-motif scores)
"""

from __future__ import annotations

import argparse
import logging
import pathlib
import time

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--abundances", nargs="+", type=float, default=[0.005, 0.01, 0.02, 0.05])
    p.add_argument("--overlaps", nargs="+", type=float, default=[0.0, 0.2, 0.4])
    p.add_argument("--methods", nargs="+", default=["harmony", "scanorama", "scvi"])
    p.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2])
    p.add_argument("--n-cells-per-batch", type=int, default=1500)
    p.add_argument("--n-batches", type=int, default=3)
    p.add_argument("--out", default="workspace/results/synth_sweep")
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    from workspace.code.datasets.synthetic import SimConfig, simulate
    from workspace.code.real_channels import ChannelConfig, loss_rate_at_k, score_purity_and_procrustes
    from workspace.code.real_channels import ot_wrapper
    from workspace.code.real_channels.ot_channel import OTChannelConfig
    from workspace.code.real_channels.nulls import ChannelNullConfig, fisher_combine, p_values_for_reports
    from workspace.code.pilots.synth_pilot import METHODS, preprocess, over_cluster, _pca_project

    out_root = pathlib.Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    rows = []
    for seed in args.seeds:
        for abundance in args.abundances:
            for overlap in args.overlaps:
                # Adjust SimConfig: abundance is 4th (epsilon); overlap overrides marker_alpha overlap
                abundances = (
                    0.45 - abundance / 3,
                    0.30 - abundance / 3,
                    0.25 - abundance / 3,
                    abundance,
                )
                cfg = SimConfig(
                    n_cells_per_batch=args.n_cells_per_batch,
                    n_batches=args.n_batches,
                    abundances=abundances,
                    seed=seed,
                )
                # patch overlap at module level
                import workspace.code.datasets.synthetic as _s
                _orig_simulate = _s.simulate

                def simulate_with_overlap(cfg=cfg, overlap=overlap):
                    # This is a cleaner hack: monkey-patch per-call
                    # But we re-implement abundance-adjusted simulate() here to control overlap
                    import anndata as ad
                    import pandas as pd
                    import scipy.sparse as sp

                    rng = np.random.default_rng(cfg.seed)
                    n_total = cfg.n_cells_per_batch * cfg.n_batches
                    labels = []
                    batches = []
                    for b in range(cfg.n_batches):
                        draw = rng.choice(cfg.n_cell_types, size=cfg.n_cells_per_batch, p=cfg.abundances)
                        labels.append(draw)
                        batches.append(np.full(cfg.n_cells_per_batch, b, dtype=np.int32))
                    labels = np.concatenate(labels)
                    batches = np.concatenate(batches)

                    markers_per_type = 100
                    base_mu = np.full((n_total, cfg.n_genes), 0.5, dtype=np.float32)
                    for t in range(cfg.n_cell_types):
                        mask = labels == t
                        base_mu[np.ix_(mask, np.arange(t * markers_per_type, (t + 1) * markers_per_type))] += cfg.marker_strength
                        if t == cfg.n_cell_types - 1:
                            base_mu[np.ix_(mask, np.arange(0, markers_per_type))] += overlap * cfg.marker_strength
                            base_mu[np.ix_(mask, np.arange(t * markers_per_type, (t + 1) * markers_per_type))] -= 0.5 * cfg.marker_strength

                    housekeeping = np.arange(cfg.n_cell_types * markers_per_type, cfg.n_genes)
                    batch_shift = rng.normal(0, cfg.batch_shift_sd, size=(cfg.n_batches, housekeeping.size)).astype(np.float32)
                    for b in range(cfg.n_batches):
                        mask = batches == b
                        base_mu[np.ix_(mask, housekeeping)] += batch_shift[b]

                    X = rng.poisson(np.exp(base_mu)).astype(np.float32)

                    obs = pd.DataFrame({
                        "batch": pd.Categorical([f"batch_{b}" for b in batches]),
                        "cell_type_public": pd.Categorical([cfg.labels[t] for t in labels]),
                        "__hidden_label__": [cfg.labels[t] for t in labels],
                        "heldout_rare": (labels == (cfg.n_cell_types - 1)),
                    })
                    obs["oracle_rare_score"] = np.where(obs["heldout_rare"], 2.0, 0.0).astype(np.float32)

                    obs_names = [f"cell_{i:06d}" for i in range(n_total)]
                    var_names = [f"g_{j:04d}" for j in range(cfg.n_genes)]
                    a = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var={"gene_ids": var_names})
                    a.obs_names = pd.Index(obs_names)
                    a.var_names = pd.Index(var_names)
                    new = a.obs["cell_type_public"].astype(object).copy()
                    new[a.obs["heldout_rare"]] = np.nan
                    a.obs["cell_type_public"] = new.astype("category")
                    return a

                a_raw = simulate_with_overlap()

                try:
                    a_pre = preprocess(a_raw, batch_key="batch", seed=seed)
                except Exception as exc:
                    logger.warning("preprocess failed for (π=%.3f, α=%.2f, s=%d): %s", abundance, overlap, seed, exc)
                    continue

                # Skip if epsilon cells got eliminated during HVG (if abundance<min)
                n_rare_cells = int(a_raw.obs["heldout_rare"].sum())
                if n_rare_cells < 3:
                    logger.warning("skipping — only %d rare cells for π=%.4f", n_rare_cells, abundance)
                    continue

                candidates = over_cluster(a_pre, seed=seed, resolution=3.0)
                hidden = a_raw.obs["__hidden_label__"].astype(str).to_numpy()
                rare_ids = set()
                is_rare = np.asarray([x == "epsilon" for x in hidden])
                for cid in np.unique(candidates):
                    m = candidates == cid
                    if m.sum() == 0:
                        continue
                    if is_rare[m].sum() >= 3 and is_rare[m].sum() / m.sum() > 0.5:
                        rare_ids.add(int(cid))

                batch_int = pd.Categorical(a_pre.obs["batch"]).codes.astype(np.int32)

                for method in args.methods:
                    if method not in METHODS:
                        continue
                    t0 = time.perf_counter()
                    try:
                        emb_int, method_meta = METHODS[method](a_pre, batch_key="batch", seed=seed)
                    except Exception as exc:
                        logger.warning("%s failed on π=%.3f α=%.2f s=%d: %s", method, abundance, overlap, seed, exc)
                        continue
                    elapsed = time.perf_counter() - t0

                    d = emb_int.shape[1]
                    pre_full = a_pre.obsm["X_uncorrected_pca"].astype(np.float32)
                    pre_trunc = pre_full[:, :d] if pre_full.shape[1] >= d else _pca_project(pre_full, d, 0)

                    cfg_ch = ChannelConfig(k=15, n_bootstrap=20, rare_abundance_threshold=max(0.10, abundance * 5), seed=seed)
                    reports = score_purity_and_procrustes(
                        emb_pre=pre_trunc,
                        emb_post=emb_int.astype(np.float32),
                        candidate_labels=candidates,
                        batch_labels=batch_int,
                        cfg=cfg_ch,
                    )
                    # OT (2-sided wrapper)
                    ot_cfg = OTChannelConfig(n_permutations=100, device="cpu", k_neighbors=15, seed=seed)
                    ot_res = ot_wrapper.score_two_sided(
                        emb_pre=pre_trunc,
                        emb_post=emb_int.astype(np.float32),
                        candidate_labels=candidates,
                        batch_labels=batch_int,
                        cfg=ot_cfg,
                    )
                    ot_by_id = {int(cid): p for cid, p in zip(ot_res["candidate_ids"], ot_res["p_values"])}

                    for r in reports:
                        cid = r.candidate_id
                        is_rare_truth = cid in rare_ids
                        rows.append({
                            "abundance": abundance,
                            "overlap": overlap,
                            "seed": seed,
                            "method": method,
                            "candidate_id": cid,
                            "candidate_abundance": r.abundance,
                            "ground_truth_label": is_rare_truth,
                            "channel_ot_p_value": ot_by_id.get(cid, float("nan")),
                            "channel_mknn_purity_drop": max(0.0, r.purity_pre - r.purity_post),
                            "channel_proc_displacement": r.procrustes_displacement,
                            "wall_clock_s": elapsed,
                        })

    df = pd.DataFrame(rows)
    df.to_parquet(out_root / "sweep_detections.parquet")
    logger.info("wrote sweep_detections.parquet (%d rows)", len(df))

    # Summary: per-(π, α, method) mean OT p-value on rare vs non-rare
    summary = df.groupby(["abundance", "overlap", "method", "ground_truth_label"]).agg(
        mean_p=("channel_ot_p_value", "mean"),
        min_p=("channel_ot_p_value", "min"),
        n=("channel_ot_p_value", "size"),
    ).reset_index()
    summary.to_csv(out_root / "sweep_summary.csv", index=False)
    logger.info("wrote sweep_summary.csv")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
