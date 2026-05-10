"""
Emit panel_A.csv + panel_B.csv for Math's ED scatter panels.

panel_A.csv: Theorem 1 rare-motif validation (n, π, Δ_σ, pred_TPR, emp_TPR).
panel_B.csv: Theorem S1-overcorrection validation (κ, Δ_σ, pred_τ², emp_detection).
"""

from __future__ import annotations

import glob
import pathlib

import numpy as np
import pandas as pd

OUT = pathlib.Path("/mnt/d-1274477442621830-m5ObBqn4/eacn3-main/examples/computational_biology/workspace/results")


def pred_TPR_theorem1(n: int, pi: float, delta_sigma: float, floor: float = 9.0, alpha: float = 0.05) -> float:
    """
    Theorem 1 predicted TPR: Φ(Z_α − √(n·π²·Δ²/(c1·log(1/αβ)·Λ))).
    We use the simplified form Φ(sqrt(signal/floor) − sqrt(1)) clipped to [0,1].
    """
    signal = n * pi * pi * (delta_sigma ** 2)
    if signal <= 0:
        return 0.0
    from scipy.stats import norm
    z_signal = np.sqrt(signal / floor)
    return float(np.clip(norm.cdf(z_signal - norm.ppf(1 - alpha)), 0.0, 1.0))


def emit_panel_A():
    """Rare-motif panel: for each (dataset, method, rare_label) compute pred vs emp TPR."""
    rows = []

    # ---------- Retina runs (real data) ----------
    # Shekhar 2016 retina via scvi.data.retina: n=19829, 2 batches.
    # Δ_σ ≈ 0.5σ per Math's envelope measurement.
    retina_abundances = {
        "retina_bc89_holdout": 0.011,
        "retina_bc4_holdout": 0.015,
        "retina_bc2_holdout": 0.021,
        "retina_bc1a_holdout": 0.042,
    }
    for holdout, pi in retina_abundances.items():
        for method in ["harmony", "scanorama", "scvi"]:
            p = f"/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/detections/real/{method}_{holdout}.parquet"
            try:
                df = pd.read_parquet(p)
            except Exception:
                continue
            rare = df[df["is_rare_truth"].astype(bool)]
            if len(rare) == 0:
                emp_tpr = 0.0
            else:
                # Emp TPR = fraction of rare-truth motifs flagged as lost (ot_p < 0.05 OR loss_prob > 0.5).
                flagged = (rare["channel_ot_p_value"] < 0.05) | (rare["loss_probability"] > 0.5)
                emp_tpr = float(flagged.mean())

            rows.append({
                "dataset": "retina",
                "method": method,
                "motif": holdout.split("_")[1],  # bc89 / bc4 / bc2 / bc1a
                "n": 19829,
                "pi": pi,
                "delta_sigma": 0.5,
                "pred_TPR": pred_TPR_theorem1(19829, pi, 0.5),
                "emp_TPR": emp_tpr,
                "color": {"harmony": "#5c8dd6", "scanorama": "#bd4848", "scvi": "#7fbd48"}.get(method, "#888"),
            })

    # ---------- Synth pancreas runs ----------
    # n=6000, epsilon at 5%, Δ_σ ≈ 1σ (measured)
    for method in ["harmony", "scanorama", "scvi"]:
        p = f"/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/detections/real/{method}_pancreas_synth.parquet"
        try:
            df = pd.read_parquet(p)
        except Exception:
            continue
        if "is_rare_truth" not in df.columns:
            # add-retroactively
            df["is_rare_truth"] = df["ground_truth_label"].astype(bool)
        rare = df[df["is_rare_truth"].astype(bool)]
        emp_tpr = float(((rare["channel_ot_p_value"] < 0.05) | (rare["loss_probability"] > 0.5)).mean()) if len(rare) else 0.0
        rows.append({
            "dataset": "pancreas_synth",
            "method": method,
            "motif": "epsilon",
            "n": 6000,
            "pi": 0.05,
            "delta_sigma": 1.0,
            "pred_TPR": pred_TPR_theorem1(6000, 0.05, 1.0),
            "emp_TPR": emp_tpr,
            "color": {"harmony": "#5c8dd6", "scanorama": "#bd4848", "scvi": "#7fbd48"}.get(method, "#888"),
        })

    df = pd.DataFrame(rows)
    out = OUT / "panel_A.csv"
    df.to_csv(out, index=False)
    print(f"Wrote panel_A.csv: {len(df)} rows")
    print(df.to_string(index=False))
    return df


def emit_panel_B():
    """Overcorrection panel: for each (dataset, method, OC_motif) compute κ, Δ_σ, emp_detection."""
    rows = []
    # Retina Scanorama/scVI overcorrection candidates (motifs 5, 11, 16, 21)
    for holdout in ["retina_bc89_holdout", "retina_bc4_holdout", "retina_bc2_holdout", "retina_bc1a_holdout"]:
        for method in ["scanorama", "scvi"]:
            p = f"/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/detections/real/{method}_{holdout}.parquet"
            try:
                df = pd.read_parquet(p)
            except Exception:
                continue
            oc = df[df.get("is_overcorrection_candidate", False).astype(bool)]
            for _, r in oc.iterrows():
                # κ = fraction of motif cells moved = 1 - bootstrap_fraction_stable approximation
                kappa = float(1.0 - r.get("channel_boot_stable", 0.5))
                # Δ_σ = Procrustes displacement normalized (rough heuristic)
                d_sigma = float(r["channel_proc_displacement"] * 100.0)  # rescale from normalized-frame
                # pred τ² = κ·Δ_σ²
                pred_tau_sq = kappa * (d_sigma ** 2)
                # emp_detection: 1 if ot_p < 0.05
                emp_detection = 1.0 if (r["channel_ot_p_value"] < 0.05) else 0.0
                rows.append({
                    "dataset": "retina",
                    "method": method,
                    "motif": f"motif_{int(r['candidate_id'])}_{holdout.split('_')[1]}",
                    "kappa": kappa,
                    "delta_sigma": d_sigma,
                    "pred_tau_sq": pred_tau_sq,
                    "emp_detection": emp_detection,
                })

    df = pd.DataFrame(rows)
    out = OUT / "panel_B.csv"
    df.to_csv(out, index=False)
    print(f"Wrote panel_B.csv: {len(df)} rows")
    print(df.head(15).to_string(index=False))
    return df


if __name__ == "__main__":
    emit_panel_A()
    print()
    emit_panel_B()
