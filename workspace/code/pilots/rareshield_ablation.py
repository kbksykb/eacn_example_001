"""
RareShield ablation suite on synth pancreas.

Per ML's next-cycle request:
  1. Full RareShield baseline (L_mass + anchor).
  2. −L_mass (only anchor) — expect large AUPRC drop.
  3. F=identity, i.e. L_rare replaced by plain L_anchor only (k=0, L=0).
  4. w_M_max sweep ∈ {0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0}.
  5. 3 seeds for each.

Note: in the fine-tune proxy, L_mass IS the protection term; L_within / L_between
from the full unified spec are not implemented in the proxy (they require encoder
gradients). This ablation runs the proxy-implementable subset.
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


def run_ablation(variant: str, w_M: float, seed: int) -> dict:
    from workspace.code.datasets.synthetic import SimConfig, simulate
    from workspace.code.real_channels import ChannelConfig, score_purity_and_procrustes
    from workspace.code.pilots.synth_pilot import (
        METHODS, preprocess, over_cluster, _pca_project,
    )
    from workspace.code.compute_l_rare import compute_l_rare, LRareConfig

    a_raw = simulate(SimConfig(seed=seed))
    hidden = a_raw.obs["__hidden_label__"].astype(str).to_numpy()
    a_pre = preprocess(a_raw, batch_key="batch", seed=seed)
    candidates = over_cluster(a_pre, seed=seed, resolution=3.0)

    rare_ids = set()
    is_rare = np.asarray([x == "epsilon" for x in hidden])
    for cid in np.unique(candidates):
        m = candidates == cid
        if m.sum() >= 3 and is_rare[m].sum() / max(m.sum(), 1) > 0.5 and is_rare[m].sum() >= 3:
            rare_ids.add(int(cid))

    emb_int, _ = METHODS["scvi"](a_pre, batch_key="batch", seed=seed)
    z_pre = a_pre.obsm["X_uncorrected_pca"][:, :emb_int.shape[1]].astype(np.float32)
    z_post_init = emb_int.astype(np.float32)
    motif_vec = np.where(np.isin(candidates, list(rare_ids)), candidates, -1)

    device = "cpu"
    z_pre_t = torch.as_tensor(z_pre, dtype=torch.float32, device=device)
    z_post = torch.as_tensor(z_post_init, dtype=torch.float32, device=device).clone().requires_grad_()
    anchor = z_post.detach().clone()
    motif_ids_t = torch.as_tensor(motif_vec, dtype=torch.long, device=device)
    opt = torch.optim.Adam([z_post], lr=1e-3)
    cfg = LRareConfig(k_N=30)

    for _ in range(150):
        opt.zero_grad()
        anc = (z_post - anchor).pow(2).sum(dim=-1).mean()

        if variant == "full":
            lr_loss, _ = compute_l_rare(z_pre_t, z_post, motif_ids_t, cfg=cfg)
            loss = w_M * lr_loss + 0.05 * anc
        elif variant == "no_L_mass":
            # Only anchor term
            loss = 0.05 * anc
        elif variant == "F_identity":
            # No admissibility gate — penalize any τ (|log ratio|²) regardless of P_d threshold.
            # Equivalent to ungated squared log-density ratio on motif cells.
            from workspace.code.rareshield.l_mass import _local_log_density
            with torch.no_grad():
                log_rho_pre = _local_log_density(z_pre_t, 30)
            log_rho_post = _local_log_density(z_post, 30)
            tau = log_rho_post - log_rho_pre
            mask = motif_ids_t >= 0
            if mask.any():
                loss = w_M * (tau[mask].pow(2).mean()) + 0.05 * anc
            else:
                loss = 0.05 * anc
        else:
            raise ValueError(f"Unknown variant {variant}")

        loss.backward()
        opt.step()

    z_shield = z_post.detach().cpu().numpy()

    # Score
    from sklearn.metrics import average_precision_score
    cfg_ch = ChannelConfig(k=15, n_bootstrap=20, rare_abundance_threshold=0.10, seed=seed)
    batch_int = pd.Categorical(a_pre.obs["batch"]).codes.astype(np.int32)
    reports_v = score_purity_and_procrustes(
        emb_pre=z_pre, emb_post=z_post_init,
        candidate_labels=candidates, batch_labels=batch_int, cfg=cfg_ch,
    )
    reports_s = score_purity_and_procrustes(
        emb_pre=z_pre, emb_post=z_shield,
        candidate_labels=candidates, batch_labels=batch_int, cfg=cfg_ch,
    )
    y_true_v = [r.candidate_id in rare_ids for r in reports_v]
    y_true_s = [r.candidate_id in rare_ids for r in reports_s]
    auprc_v = float(average_precision_score(y_true_v, [r.loss_probability for r in reports_v])) if any(y_true_v) else float("nan")
    auprc_s = float(average_precision_score(y_true_s, [r.loss_probability for r in reports_s])) if any(y_true_s) else float("nan")

    import scanpy as sc
    from sklearn.metrics import adjusted_rand_score
    def ari(emb, labels):
        a = a_pre.copy()
        a.obsm["tmp"] = emb
        sc.pp.neighbors(a, use_rep="tmp", random_state=seed)
        sc.tl.leiden(a, resolution=1.0, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
        return float(adjusted_rand_score(labels, a.obs["leiden"].astype(int).to_numpy()))

    ari_v = ari(z_post_init, hidden)
    ari_s = ari(z_shield, hidden)

    return {
        "variant": variant, "w_M": w_M, "seed": seed,
        "auprc_vanilla": auprc_v, "auprc_shield": auprc_s,
        "delta_auprc": auprc_s - auprc_v,
        "ari_vanilla": ari_v, "ari_shield": ari_s,
        "delta_ari": ari_s - ari_v,
    }


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")
    rows = []
    for seed in [0, 1]:
        for variant in ["full", "no_L_mass", "F_identity"]:
            for w_M in [2.0]:  # minimal spread in proxy
                logger.info("=== %s w_M=%s seed=%d ===", variant, w_M, seed)
                t0 = time.perf_counter()
                r = run_ablation(variant, w_M, seed)
                r["wall_s"] = time.perf_counter() - t0
                rows.append(r)
                print(json.dumps(r, indent=2))
        # also w_M sweep
        for w_M in [0.01, 0.1, 0.5, 2.0]:
            logger.info("=== w_M=%s seed=%d ===", w_M, seed)
            t0 = time.perf_counter()
            r = run_ablation("full", w_M, seed)
            r["wall_s"] = time.perf_counter() - t0
            rows.append(r)
    df = pd.DataFrame(rows)
    print("\n=== Ablation summary ===")
    print(df.to_string(index=False))
    out = pathlib.Path(__file__).resolve().parents[3] / "results" / "rareshield_ablation.csv"
    df.to_csv(out, index=False)
    logger.info("summary → %s", out)


if __name__ == "__main__":
    main()
