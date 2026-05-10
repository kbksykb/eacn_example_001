"""Inspect BBKNN blind results vs pre-registered prediction.

Reads the .parquet and .npy written by bbknn_prospective_run.py, compares to
the registered prediction:
    κ_med ≤ 2.5 on both cohorts
    Cheng-PAAD p_bh ≥ 0.5
    Cheng6 p_bh < 0.01

Writes `workspace/results/bbknn_prospective_result.md` with pass/fail per
criterion and the full numbers.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

SHARED = Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared")


def main() -> None:
    rows = []
    for name in ["cheng_paad_lamp3_holdout", "cheng6_lamp3_holdout"]:
        det_path = SHARED / "detections" / "real" / f"bbknn_{name}.parquet"
        k_path = SHARED / "detections" / "real" / f"bbknn_{name}_kappa.npy"
        df = pd.read_parquet(det_path)
        kappas = np.load(k_path)

        rare = df[df["ground_truth_label"]]
        non_rare = df[~df["ground_truth_label"]]
        if len(rare) == 0:
            rare_min_p_raw = 1.0
        else:
            rare_min_p_raw = float(rare["channel_ot_p_value"].min())
        # BH per run
        from scipy.stats import false_discovery_control
        raw_p_all = df["channel_ot_p_value"].to_numpy()
        bh_p_all = false_discovery_control(raw_p_all, method="bh")
        id_to_bh = dict(zip(df["candidate_id"].to_numpy(), bh_p_all))
        rare_ids = rare["candidate_id"].tolist()
        rare_min_p_bh = min([id_to_bh[c] for c in rare_ids]) if rare_ids else 1.0
        non_rare_flag_bh = float(sum(1 for p in [id_to_bh[c] for c in non_rare["candidate_id"]] if p < 0.05) / max(len(non_rare), 1))

        k_med = float(np.median(kappas)) if len(kappas) else 0.0
        k_q1 = float(np.percentile(kappas, 25)) if len(kappas) else 0.0
        k_q3 = float(np.percentile(kappas, 75)) if len(kappas) else 0.0

        rows.append(dict(
            dataset=name,
            n_candidates=len(df),
            n_rare_motifs=len(rare_ids),
            rare_min_p_raw=rare_min_p_raw,
            rare_min_p_bh=rare_min_p_bh,
            non_rare_flag_bh=non_rare_flag_bh,
            kappa_median=k_med, kappa_q1=k_q1, kappa_q3=k_q3,
            n_kappa_motifs=len(kappas),
        ))

    # Evaluate against prediction
    paad = next(r for r in rows if "paad" in r["dataset"])
    cheng6 = next(r for r in rows if r["dataset"] == "cheng6_lamp3_holdout")

    # Criterion evaluation
    paad_k_ok = paad["kappa_median"] <= 2.5
    cheng6_k_ok = cheng6["kappa_median"] <= 2.5
    paad_preserve_ok = paad["rare_min_p_bh"] >= 0.5
    cheng6_fire_ok = cheng6["rare_min_p_bh"] < 0.01

    all_confirm = all([paad_k_ok, cheng6_k_ok, paad_preserve_ok, cheng6_fire_ok])

    if all_confirm:
        outcome_txt = "All four pre-registered criteria pass. BBKNN exhibits Regime A (low-κ, scale-monotonic) exactly as Math predicted. The κ three-regime mechanism is predictive, not just retrospective."
        next_txt = ("- BBKNN κ three-regime prediction holds prospectively. Strengthens the paper considerably — the mechanism is not just curve-fit to Harmony/Scanorama/scVI, it predicts 4th-integrator behavior correctly.\n"
                    "- Extend κ grid to 17 points (9 Cheng + 3 Zheng + 3 Sade + 2 BBKNN).\n"
                    "- Math to cite this in §rem:S1-three-regime-kappa.\n"
                    "- Consider 5th integrator (e.g. scANVI) as second prospective test.")
    else:
        outcome_txt = "At least one pre-registered criterion failed. See criteria table. This disconfirms either the κ bucket assignment for BBKNN specifically, or the universality of the three-regime taxonomy. Next step: document the failure mode in Supp S5."
        next_txt = ("- Document failure in Supp S5.2 alongside ℓ-unification disconfirmation.\n"
                    "- Investigate which criterion failed and mechanistic explanation.\n"
                    "- Consider refinement: does BBKNN land in a genuinely new regime, or is the κ bucket threshold wrong?")

    doc = f"""# BBKNN prospective prediction — result

**Date:** 2026-05-11
**Pre-registration:** `workspace/results/bbknn_prospective_prediction.md` (Math's formal prediction, committed before data observation)
**Status:** BLIND RUN COMPLETE. Results now inspected.

## Result numbers

| Dataset | κ_median | IQR | rare_min_p_raw | rare_min_p_bh | non_rare_flag_bh | #rare_motifs |
|---------|---------:|----:|---------------:|--------------:|-----------------:|------------:|
| Cheng PAAD 2.8k | {paad['kappa_median']:.3f} | [{paad['kappa_q1']:.3f}, {paad['kappa_q3']:.3f}] | {paad['rare_min_p_raw']:.4f} | {paad['rare_min_p_bh']:.4f} | {paad['non_rare_flag_bh']:.3f} | {paad['n_rare_motifs']} |
| Cheng6 49k | {cheng6['kappa_median']:.3f} | [{cheng6['kappa_q1']:.3f}, {cheng6['kappa_q3']:.3f}] | {cheng6['rare_min_p_raw']:.4f} | {cheng6['rare_min_p_bh']:.4f} | {cheng6['non_rare_flag_bh']:.3f} | {cheng6['n_rare_motifs']} |

## Prediction evaluation

| Criterion | Prediction | Observation | Pass? |
|-----------|------------|-------------|:-----:|
| κ_med Cheng PAAD ≤ 2.5 | ≤ 2.5 | {paad['kappa_median']:.3f} | {'✓' if paad_k_ok else '✗'} |
| κ_med Cheng6 ≤ 2.5 | ≤ 2.5 | {cheng6['kappa_median']:.3f} | {'✓' if cheng6_k_ok else '✗'} |
| Cheng PAAD p_bh ≥ 0.5 (preserve) | ≥ 0.5 | {paad['rare_min_p_bh']:.4f} | {'✓' if paad_preserve_ok else '✗'} |
| Cheng6 p_bh < 0.01 (fire) | < 0.01 | {cheng6['rare_min_p_bh']:.4f} | {'✓' if cheng6_fire_ok else '✗'} |

## Outcome

**{'FULL CONFIRMATION' if all_confirm else 'DISCONFIRMATION (see criteria table)'}**

{outcome_txt}

## Register-compliance notes

- Prediction was registered BEFORE run (commit da01cac).
- Runner was blind until this evaluation file was written.
- Disconfirmation conditions were pre-specified.
- Result is reported at face value; no re-fitting of criteria post-hoc.

## Next steps

{next_txt}
"""
    out = Path("workspace/results/bbknn_prospective_result.md")
    out.write_text(doc)
    print(f"wrote {out}")
    print(f"CONFIRMED = {all_confirm}")


if __name__ == "__main__":
    main()
