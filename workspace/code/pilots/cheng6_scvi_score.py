"""Enrich Cheng6 scVI integration with motif_id_pre, __hidden_label__, X_uncorrected_pca
copied from the Harmony integration (same obs_names order) so the REAL scoring
pipeline (panel_A_tight.py) can run against it.

Then score via rerun_ot_on_h5ad at B=500, BH-FDR within dataset, and emit
the parquet + Panel-A-row inline.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

SHARED = Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared")


def enrich() -> Path:
    src = ad.read_h5ad(SHARED / "integrations" / "scvi" / "cheng6_lamp3_holdout.h5ad")
    ref = ad.read_h5ad(SHARED / "integrations" / "harmony" / "cheng6_lamp3_holdout.h5ad")
    assert list(src.obs_names) == list(ref.obs_names), "obs_names order mismatch between scvi and harmony outputs"

    src.obsm["X_integrated"] = src.obsm["X_scvi"].astype(np.float32)
    src.obsm["X_uncorrected_pca"] = ref.obsm["X_uncorrected_pca"].astype(np.float32)

    for key in ["motif_id_pre", "__hidden_label__", "heldout_rare", "cell_type_public",
                "leiden", "oracle_rare_score", "MajorCluster", "labels"]:
        if key in ref.obs and key not in src.obs:
            src.obs[key] = ref.obs[key].values

    out = SHARED / "integrations" / "scvi" / "cheng6_lamp3_holdout.h5ad"
    src.write_h5ad(out)
    print(f"[enrich] wrote enriched {out}  obs_cols has motif_id_pre={'motif_id_pre' in src.obs}")
    return out


def score(h5ad: Path, n_perm: int = 500, seed: int = 0):
    sys.path.insert(0, "/mnt/d-1274477442621830-m5ObBqn4/eacn3-main/examples/computational_biology")
    from workspace.code.real_channels import ot_wrapper
    from workspace.code.real_channels.ot_channel import OTChannelConfig

    a = ad.read_h5ad(h5ad)
    z_pre = a.obsm["X_uncorrected_pca"].astype(np.float32)
    z_post = a.obsm["X_integrated"].astype(np.float32)
    d = z_post.shape[1]
    pre_trunc = z_pre[:, :d]
    motif_ids = a.obs["motif_id_pre"].astype(int).to_numpy()
    batch_int = pd.Categorical(a.obs["batch"]).codes.astype(np.int32)

    print(f"[score] OT B={n_perm} on {h5ad.name}  n={a.n_obs}  d_post={d}  motifs={len(np.unique(motif_ids))}")
    cfg = OTChannelConfig(n_permutations=n_perm, device="cpu", k_neighbors=15, seed=seed)
    r = ot_wrapper.score_two_sided(pre_trunc, z_post, motif_ids, batch_int, cfg)
    cand_ids = r["candidate_ids"].tolist()
    raw_p = r["p_values"].tolist()
    stat = r["stat"].tolist()
    signed = r["signed_stat_mean"].tolist()

    bh_p = list(false_discovery_control(np.asarray(raw_p), method="bh"))

    hidden = a.obs["__hidden_label__"].astype(str).to_numpy()
    is_lamp3 = np.array([("LAMP3" in x) or ("cDC3" in x) or (x == "M05_cDC3_LAMP3") for x in hidden])

    rare_ids: set[int] = set()
    for cid in cand_ids:
        mask = motif_ids == cid
        if mask.sum() == 0:
            continue
        hit_frac = is_lamp3[mask].sum() / mask.sum()
        if is_lamp3[mask].sum() >= 3 and hit_frac > 0.5:
            rare_ids.add(int(cid))
    print(f"[score] rare (LAMP3+ majority motif) candidate_ids: {sorted(rare_ids)}  n={len(rare_ids)}")

    rows = []
    for cid, st, rp, bp, ss in zip(cand_ids, stat, raw_p, bh_p, signed):
        cid = int(cid)
        is_rare = cid in rare_ids
        abundance = float((motif_ids == cid).mean())
        rows.append(dict(
            method="scvi", dataset="cheng6_lamp3_holdout", candidate_id=cid,
            abundance=abundance,
            channel_mknn_purity_pre=np.nan, channel_mknn_purity_post=np.nan,
            channel_proc_displacement=np.nan, channel_boot_stable=np.nan,
            channel_ot_stat=float(st), channel_ot_p_value=float(rp), channel_ot_signed=float(ss),
            loss_probability=float(np.clip(1.0 - bp, 0.0, 1.0)),
            ground_truth_label=is_rare, seed=seed, n_cells=int(a.n_obs),
        ))
    df = pd.DataFrame(rows)
    out = SHARED / "detections" / "real" / "scvi_cheng6_lamp3_holdout.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False)
    print(f"[score] wrote {out}  rows={len(df)}")

    rare_min_raw = df.loc[df.ground_truth_label, "channel_ot_p_value"].min() if df.ground_truth_label.any() else 1.0
    rare_min_bh = df.loc[df.ground_truth_label, "loss_probability"].apply(lambda x: 1 - x).min() if df.ground_truth_label.any() else 1.0
    # Actually bh p-values: reconstruct
    bh_p_map = dict(zip(cand_ids, bh_p))
    rare_min_bh = min([bh_p_map[c] for c in rare_ids]) if rare_ids else 1.0
    nr_flag_raw = ((~df.ground_truth_label) & (df.channel_ot_p_value < 0.05)).sum() / max((~df.ground_truth_label).sum(), 1)
    nr_flag_bh_ids = [c for c in cand_ids if c not in rare_ids and bh_p_map[c] < 0.05]
    nr_flag_bh = len(nr_flag_bh_ids) / max(len(cand_ids) - len(rare_ids), 1)

    print(f"[score] rare_min_p_raw={rare_min_raw:.4f}  rare_min_p_bh={rare_min_bh:.4f}  non_rare_flag_bh={nr_flag_bh:.3f}")
    return dict(rare_min_p_raw=rare_min_raw, rare_min_p_bh=rare_min_bh,
                non_rare_flag_raw=nr_flag_raw, non_rare_flag_bh=nr_flag_bh,
                n_total=len(cand_ids), n_rare=len(rare_ids))


def append_panel_A(stats):
    import csv
    path = Path("workspace/results/panel_A_tight_v2.csv")
    # Compute matches per schema from tail of panel_A_tight_v2:
    # cheng6_pancancer,cheng6_lamp3,harmony,49271,0.0115,0.8,0.001996...,0.004109...,0.4706,0.4706,0.0,0.0,35,1,1.0
    flagged_strict = 1.0 if (stats["rare_min_p_bh"] < 0.05) else 0.0
    row = [
        "cheng6_pancancer", "cheng6_lamp3", "scvi",
        49271, 0.0115, 0.8,
        round(stats["rare_min_p_raw"], 14), round(stats["rare_min_p_bh"], 14),
        round(stats["non_rare_flag_raw"], 4), round(stats["non_rare_flag_bh"], 4),
        0.0, 0.0,
        stats["n_total"] - stats["n_rare"], stats["n_rare"],
        flagged_strict,
    ]
    with path.open("a", newline="") as f:
        csv.writer(f).writerow(row)
    print(f"[panel_A] appended row: {row}")


def main() -> None:
    h5ad = enrich()
    stats = score(h5ad, n_perm=500)
    append_panel_A(stats)


if __name__ == "__main__":
    main()
