"""
Recompute Panel A with tight TPR: B=500 permutation null + BH-FDR per-holdout,
rare-motif flagged AND non-rare specificity > 0.95.

Emits:
    workspace/results/panel_A_tight_v2.csv with columns:
        dataset, holdout, method, n, pi, delta_sigma,
        raw_tight_TPR (no FDR), bh_tight_TPR (BH-FDR applied),
        rare_min_p_raw, rare_min_p_bh, non_rare_flag_rate_raw, non_rare_flag_rate_bh
"""

from __future__ import annotations

import logging
import pathlib

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

logger = logging.getLogger(__name__)


def rerun_ot_on_h5ad(h5ad_path: str, n_perm: int = 500, seed: int = 0):
    """Rerun OT channel with B=500 on a specific integration H5AD."""
    import anndata as ad

    a = ad.read_h5ad(h5ad_path)
    from workspace.code.real_channels import ot_wrapper
    from workspace.code.real_channels.ot_channel import OTChannelConfig

    z_pre = a.obsm["X_uncorrected_pca"].astype(np.float32)
    z_post = a.obsm["X_integrated"].astype(np.float32)
    d = z_post.shape[1]
    pre_trunc = z_pre[:, :d]
    motif_ids = a.obs["motif_id_pre"].astype(int).to_numpy()
    batch_int = pd.Categorical(a.obs["batch"]).codes.astype(np.int32)

    cfg = OTChannelConfig(n_permutations=n_perm, device="cpu", k_neighbors=15, seed=seed)
    r = ot_wrapper.score_two_sided(pre_trunc, z_post, motif_ids, batch_int, cfg)
    return dict(zip(r["candidate_ids"].tolist(), r["p_values"].tolist())), a


def bh_adjust(p_values: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR adjustment."""
    p = np.asarray(p_values)
    return list(false_discovery_control(p, method="bh"))


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    configs = []
    # Retina holdouts × 3 methods
    for holdout in ["retina_bc89_holdout", "retina_bc4_holdout", "retina_bc2_holdout", "retina_bc1a_holdout"]:
        abu = {"bc89": 0.011, "bc4": 0.015, "bc2": 0.021, "bc1a": 0.042}[holdout.split("_")[1]]
        for method in ["harmony", "scanorama", "scvi"]:
            configs.append(("retina", holdout, method, 19829, abu, 0.5, "BC8_9" if holdout=="retina_bc89_holdout" else holdout.split("_")[1].upper()))
    # Synth pancreas
    for method in ["harmony", "scanorama", "scvi"]:
        configs.append(("pancreas_synth", "pancreas_synth", method, 6000, 0.05, 1.0, "epsilon"))
    # Cheng5 pan-cancer LAMP3+
    for method in ["harmony", "scanorama", "scvi"]:
        configs.append(("cheng5_pancancer", "cheng5_lamp3_holdout", method, 20341, 0.015, 0.8, "LAMP3_mregDC"))

    rows = []
    for dataset, holdout, method, n, pi, delta_sigma, rare_label in configs:
        h5ad = f"/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/integrations/{method}/{holdout}.h5ad"
        if not pathlib.Path(h5ad).exists():
            logger.warning("skip — %s not found", h5ad)
            continue
        logger.info("rerunning OT B=500 on %s/%s", method, holdout)
        ot_by_id, a = rerun_ot_on_h5ad(h5ad)

        # Figure out rare-truth motif ids on the fly
        hidden = a.obs.get("__hidden_label__", pd.Series(["?"] * a.n_obs)).astype(str).to_numpy()
        candidates = a.obs["motif_id_pre"].astype(int).to_numpy()
        is_rare = np.array([str(x).upper().replace("_", "").replace("-","") == rare_label.upper().replace("_", "").replace("-","") for x in hidden])
        rare_ids = set()
        for cid in np.unique(candidates):
            m = candidates == cid
            if m.sum() == 0: continue
            if is_rare[m].sum() >= 3 and is_rare[m].sum() / m.sum() > 0.5:
                rare_ids.add(int(cid))

        # p-values for all motifs
        motifs_sorted = sorted(ot_by_id)
        p_raw = [ot_by_id[cid] for cid in motifs_sorted]
        p_bh = bh_adjust(p_raw)

        rare_p_raw = [p for cid, p in zip(motifs_sorted, p_raw) if cid in rare_ids]
        rare_p_bh = [p for cid, p in zip(motifs_sorted, p_bh) if cid in rare_ids]
        non_rare_p_raw = [p for cid, p in zip(motifs_sorted, p_raw) if cid not in rare_ids]
        non_rare_p_bh = [p for cid, p in zip(motifs_sorted, p_bh) if cid not in rare_ids]

        def tpr(rare_p, non_rare_p, alpha=0.05, spec_thresh=0.95):
            if not rare_p:
                return 0.0
            rare_flagged = any(p < alpha for p in rare_p)
            non_rare_flag_rate = sum(p < alpha for p in non_rare_p) / max(len(non_rare_p), 1)
            specific = (1 - non_rare_flag_rate) > spec_thresh
            return 1.0 if (rare_flagged and specific) else 0.0

        rows.append({
            "dataset": dataset,
            "holdout": holdout.replace("retina_", "").replace("_holdout", ""),
            "method": method,
            "n": n, "pi": pi, "delta_sigma": delta_sigma,
            "rare_min_p_raw": min(rare_p_raw) if rare_p_raw else float("nan"),
            "rare_min_p_bh": min(rare_p_bh) if rare_p_bh else float("nan"),
            "non_rare_flag_rate_raw": sum(p < 0.05 for p in non_rare_p_raw) / max(len(non_rare_p_raw), 1),
            "non_rare_flag_rate_bh": sum(p < 0.05 for p in non_rare_p_bh) / max(len(non_rare_p_bh), 1),
            "raw_tight_TPR": tpr(rare_p_raw, non_rare_p_raw),
            "bh_tight_TPR": tpr(rare_p_bh, non_rare_p_bh),
            "n_motifs": len(motifs_sorted),
            "n_rare_motifs": len(rare_ids),
        })

    df = pd.DataFrame(rows)
    out = pathlib.Path("/mnt/d-1274477442621830-m5ObBqn4/eacn3-main/examples/computational_biology/workspace/results/panel_A_tight_v2.csv")
    df.to_csv(out, index=False)
    print(df.to_string(index=False))
    logger.info("wrote %s", out)


if __name__ == "__main__":
    main()
