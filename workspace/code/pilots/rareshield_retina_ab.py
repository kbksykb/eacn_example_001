"""
RareShield A/B on retina integrated embeddings.

For each method in {harmony, scanorama, scvi}:
    1. Load /mnt/.../shared/integrations/<method>/retina_bc89_holdout.h5ad
       (already has obsm['X_uncorrected_pca'] pre, obsm['X_integrated'] post)
    2. Build motif_id_pre from the existing obs column.
    3. Identify rare-truth motifs (BC8_9 majority) and overcorrection-candidates
       (non-rare motifs with loss_probability > 0.5 from the detection parquet).
    4. Fine-tune obsm['X_integrated'] against Math's compute_l_rare (L_mass) with
       motif_ids pointing at BOTH rare-truth AND overcorrection-candidates.
    5. Re-score REAL channels on (pre vs post_shielded).
    6. Report:
       - AUPRC on rare-truth (preservation of rare).
       - Overcorrection-flag p-value drop on abundants (protection against
         overcorrection).
       - ARI on true retina labels (1-based unsuprevised Leiden recoveries).

This is the retina equivalent of workspace/code/pilots/rareshield_ab.py (synth).
"""

from __future__ import annotations

import argparse
import json
import logging
import pathlib
import time

import numpy as np
import pandas as pd
import torch

logger = logging.getLogger(__name__)

SHARED = pathlib.Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared")


def fine_tune(
    z_pre: np.ndarray,
    z_post_init: np.ndarray,
    motif_ids: np.ndarray,
    lambda_mass: float = 2.0,
    n_epochs: int = 100,
    lr: float = 1e-3,
    seed: int = 0,
    subsample_n: int | None = 4000,
):
    from workspace.code.compute_l_rare import compute_l_rare, LRareConfig
    torch.manual_seed(seed)
    device = "cpu"  # pairwise-dist O(n²) kills GPU at n=20k; CPU is fine for the small fine-tune
    rng = np.random.default_rng(seed)

    # If data is big, subsample — but always include all motif cells
    n_total = z_pre.shape[0]
    if subsample_n is not None and n_total > subsample_n:
        motif_idx = np.where(motif_ids >= 0)[0]
        keep_motif = motif_idx
        other_pool = np.setdiff1d(np.arange(n_total), motif_idx)
        n_other = max(subsample_n - len(motif_idx), 0)
        keep_other = rng.choice(other_pool, size=min(n_other, len(other_pool)), replace=False)
        keep = np.sort(np.concatenate([keep_motif, keep_other]))
    else:
        keep = np.arange(n_total)

    z_pre_sub = z_pre[keep]
    z_post_sub = z_post_init[keep].copy()
    motif_sub = motif_ids[keep]

    z_pre_t = torch.as_tensor(z_pre_sub, dtype=torch.float32, device=device)
    z_post = torch.as_tensor(z_post_sub, dtype=torch.float32, device=device).requires_grad_()
    anchor = z_post.detach().clone()
    motif_ids_t = torch.as_tensor(motif_sub, dtype=torch.long, device=device)
    opt = torch.optim.Adam([z_post], lr=lr)
    cfg = LRareConfig(k_N=30)
    trace = []
    for ep in range(n_epochs):
        opt.zero_grad()
        lr_loss, _ = compute_l_rare(z_pre_t, z_post, motif_ids_t, cfg=cfg)
        anc = (z_post - anchor).pow(2).sum(dim=-1).mean()
        loss = lambda_mass * lr_loss + 0.05 * anc
        loss.backward()
        opt.step()
        if ep % 25 == 0:
            trace.append({"ep": ep, "loss": float(loss.detach().cpu()),
                         "lr": float(lr_loss.detach().cpu()),
                         "anchor": float(anc.detach().cpu())})

    # Apply the learned perturbation to the subsampled cells; non-subsampled stay vanilla.
    z_out = z_post_init.copy()
    z_out[keep] = z_post.detach().cpu().numpy()
    return z_out, trace


def score_and_ari(adata_pre, z_pre, z_post, motif_candidates, batch_int, true_labels):
    from workspace.code.real_channels import ChannelConfig, score_purity_and_procrustes
    from workspace.code.real_channels import ot_wrapper
    from workspace.code.real_channels.ot_channel import OTChannelConfig
    from sklearn.metrics import adjusted_rand_score
    import scanpy as sc
    cfg = ChannelConfig(k=15, n_bootstrap=20, rare_abundance_threshold=0.10, seed=0)
    reports = score_purity_and_procrustes(
        emb_pre=z_pre[:, :z_post.shape[1]].astype(np.float32),
        emb_post=z_post.astype(np.float32),
        candidate_labels=motif_candidates,
        batch_labels=batch_int,
        cfg=cfg,
    )
    ot_cfg = OTChannelConfig(n_permutations=60, device="cpu", k_neighbors=15, seed=0)
    ot_res = ot_wrapper.score_two_sided(
        emb_pre=z_pre[:, :z_post.shape[1]].astype(np.float32),
        emb_post=z_post.astype(np.float32),
        candidate_labels=motif_candidates,
        batch_labels=batch_int,
        cfg=ot_cfg,
    )
    ot_by_id = dict(zip(ot_res["candidate_ids"].tolist(), ot_res["p_values"].tolist()))

    a = adata_pre.copy()
    a.obsm["__tmp"] = z_post
    sc.pp.neighbors(a, use_rep="__tmp", random_state=0)
    sc.tl.leiden(a, resolution=1.0, random_state=0, flavor="igraph", n_iterations=2, directed=False)
    ari = adjusted_rand_score(true_labels, a.obs["leiden"].astype(int).to_numpy())

    return reports, ot_by_id, float(ari)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--method", default="scanorama", choices=["harmony", "scanorama", "scvi"])
    p.add_argument("--rare-label", default="BC8_9")
    p.add_argument("--holdout-name", default="retina_bc89_holdout")
    p.add_argument("--lambda-mass", type=float, default=2.0)
    p.add_argument("--include-overcorrection-motifs", action="store_true",
                   help="Include motifs flagged as overcorrected (loss_prob > 0.5) as L_mass targets")
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")

    import anndata as ad
    h5ad_path = SHARED / "integrations" / args.method / f"{args.holdout_name}.h5ad"
    a = ad.read_h5ad(h5ad_path)
    logger.info("loaded %s shape=%s", h5ad_path, a.shape)

    z_pre = a.obsm["X_uncorrected_pca"].astype(np.float32)
    z_post = a.obsm["X_integrated"].astype(np.float32)
    candidates = a.obs["motif_id_pre"].astype(int).to_numpy()
    batch_int = pd.Categorical(a.obs["batch"]).codes.astype(np.int32)

    # Identify rare-truth motifs
    hidden = a.obs["__hidden_label__"].astype(str).to_numpy()
    is_rare = np.asarray([x == args.rare_label for x in hidden])
    rare_ids = set()
    for cid in np.unique(candidates):
        m = candidates == cid
        if m.sum() >= 3 and is_rare[m].sum() >= 3 and is_rare[m].sum() / m.sum() > 0.5:
            rare_ids.add(int(cid))

    # Optionally include overcorrection-candidates from the detection parquet
    oc_ids: set[int] = set()
    if args.include_overcorrection_motifs:
        det_path = SHARED / "detections" / "real" / f"{args.method}_{args.holdout_name}.parquet"
        if det_path.exists():
            df = pd.read_parquet(det_path)
            oc_ids = set(df.loc[(df["loss_probability"] > 0.5) & (~df["ground_truth_label"]), "candidate_id"].astype(int))
            logger.info("overcorrection-candidate motif ids: %s", sorted(oc_ids))

    target_ids = rare_ids | oc_ids
    logger.info("L_mass target motif ids: %s", sorted(target_ids))
    motif_vec = np.where(np.isin(candidates, list(target_ids)), candidates, -1)

    # Truncate pre to same d as post
    pre_trunc = z_pre[:, :z_post.shape[1]]
    # Vanilla rescoring (for baseline)
    reports_v, ot_v, ari_v = score_and_ari(a, pre_trunc, z_post, candidates, batch_int, hidden)

    # Fine-tune
    logger.info("fine-tuning %s with λ=%s…", args.method, args.lambda_mass)
    t0 = time.perf_counter()
    z_shield, trace = fine_tune(pre_trunc, z_post, motif_vec, lambda_mass=args.lambda_mass)
    t_shield = time.perf_counter() - t0
    logger.info("fine-tune done in %.1fs", t_shield)
    reports_s, ot_s, ari_s = score_and_ari(a, pre_trunc, z_shield, candidates, batch_int, hidden)

    # Build comparison table
    rows_v = [(r.candidate_id, r.candidate_id in rare_ids, r.loss_probability, ot_v.get(r.candidate_id, np.nan)) for r in reports_v]
    rows_s = [(r.candidate_id, r.candidate_id in rare_ids, r.loss_probability, ot_s.get(r.candidate_id, np.nan)) for r in reports_s]
    df_v = pd.DataFrame(rows_v, columns=["cid", "gt", "lp_v", "ot_v"])
    df_s = pd.DataFrame(rows_s, columns=["cid", "gt", "lp_s", "ot_s"])
    joined = df_v.merge(df_s, on="cid", how="outer")
    joined["gt"] = joined["gt_x"].fillna(joined["gt_y"])
    print(joined[["cid", "gt", "lp_v", "lp_s", "ot_v", "ot_s"]].to_string(index=False))

    # Overcorrection-flag p-value drop: average OT p-value on OC candidates before/after
    oc_mask = joined["cid"].isin(list(oc_ids))
    rare_mask = joined["gt"].astype(bool)
    ot_drop_on_oc = float((joined.loc[oc_mask, "ot_s"] - joined.loc[oc_mask, "ot_v"]).mean()) if oc_mask.any() else float("nan")
    ot_drop_on_rare = float((joined.loc[rare_mask, "ot_s"] - joined.loc[rare_mask, "ot_v"]).mean()) if rare_mask.any() else float("nan")
    # Loss-probability changes
    lp_drop_on_oc = float((joined.loc[oc_mask, "lp_s"] - joined.loc[oc_mask, "lp_v"]).mean()) if oc_mask.any() else float("nan")
    lp_drop_on_rare = float((joined.loc[rare_mask, "lp_s"] - joined.loc[rare_mask, "lp_v"]).mean()) if rare_mask.any() else float("nan")

    summary = {
        "method": args.method,
        "holdout_name": args.holdout_name,
        "lambda_mass": args.lambda_mass,
        "include_overcorrection_motifs": args.include_overcorrection_motifs,
        "n_rare_truth_motifs": len(rare_ids),
        "n_overcorrection_target_motifs": len(oc_ids),
        "ARI_vanilla": ari_v,
        "ARI_rareshield": ari_s,
        "delta_ARI": ari_s - ari_v,
        "mean_OT_p_delta_rare_truth": ot_drop_on_rare,
        "mean_OT_p_delta_overcorrection_candidate": ot_drop_on_oc,
        "mean_lp_delta_rare_truth": lp_drop_on_rare,
        "mean_lp_delta_overcorrection_candidate": lp_drop_on_oc,
        "fine_tune_seconds": t_shield,
        "trace": trace,
    }
    print(json.dumps({k: v for k, v in summary.items() if k != "trace"}, indent=2))

    out = pathlib.Path(__file__).resolve().parents[3] / "results" / f"rareshield_retina_{args.method}_{args.holdout_name}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2))
    logger.info("summary → %s", out)


if __name__ == "__main__":
    main()
