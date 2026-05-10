# HVG3K preprocessing effect on κ — Candidate A vs Candidate B discrimination

**Date:** 2026-05-11
**Input:** The Sade-Feldman κ disconfirmation (κ_scVI_Sade_HVG3K=1.15, LOWER than full-gene Cheng6 scVI κ=5.13 and reversing the Harmony < scVI ordering) had two candidate explanations:
- **Candidate A**: HVG3K preprocessing shrinks κ (preprocessing artifact).
- **Candidate B**: State-regime rare types (fragmented across motifs) genuinely break the κ mechanism.

**Test protocol**: Rerun Cheng6 Scanorama and scVI on HVG3K preprocessing, measure κ, compare to full-gene Cheng6 κ values. Same dataset, same rare type (LAMP3+ mregDC), only preprocessing differs.

## Result — both methods now tested

| Integrator | Preprocessing | n × genes | κ_median | IQR |
|-----------|---------------|----------:|---------:|----:|
| scVI | **Full-gene** | 49,271 × 11,104 | **5.13** | [4.24, 6.95] |
| scVI | **HVG3K** | 49,271 × 3,000 | **6.44** | [5.80, 7.70] |
| Scanorama | Full-gene | 49,271 × 11,104 | **9.13** | [7.58, 10.99] |
| Scanorama | HVG3K | 49,271 × 3,000 | **9.82** | [8.82, 11.34] |

## Candidate A disconfirmed for both methods

HVG3K preprocessing **slightly increases** κ for both scVI (5.13 → 6.44, +25%) and Scanorama (9.13 → 9.82, +8%). The direction is consistently upward, not downward. **Candidate A is therefore disconfirmed** — HVG3K preprocessing does NOT shrink κ on Cheng6 data.

Consequence: the Sade scVI κ=1.15 and Sade Scanorama HVG3K κ=3.89 (both LOWER than full-gene Cheng6 values) cannot be attributed to HVG3K preprocessing. The shrinkage must come from something about Sade-Feldman itself.

## Candidate B now has positive support

Sade-Feldman canonical TPEX is a **fragmented state** distributed across 6 motifs at 8.9-16.8% enrichment each (per workspace/results/sade_tpex_pilot.md), not a discrete cluster. Motif-centroid κ is the wrong statistic for this regime:
- When the rare cells are distributed across motifs, no single motif carries enough rare-cell mass for centroid displacement to track rare-cell redistribution faithfully.
- scVI's low κ=1.15 on Sade reflects scVI's mild regularisation preserving the OVERALL state geometry (where cells with TPEX signature are still near each other across batches), but across-motif centroid displacement is insensitive to whether the rare TPEX STATE is preserved or fragmented further.

## Scope refinement for §4.6

The κ three-regime framework is scoped to **cluster-regime rare types** (discrete populations with coherent motif assignment). It is **not well-defined for state-regime rare types** (distributed across motifs) — for those, κ measured on motif centroids does not track the rare-state-preservation signal.

This is another scope caveat to disclose alongside the adapter-specification caveat for graph-only integrators.

## What this means for the paper

- κ three-regime stratification remains confirmed on **cluster-regime full-gene Cheng + Zheng MM data** (12 points, Harmony < scVI < Scanorama ordering preserved exactly).
- κ three-regime stratification does NOT generalize to **state-regime Sade TPEX data**, and the scope limit is documented.
- Detection (REAL) still works on state-regime rare types (per sade_tpex_pilot.md; Harmony + Scanorama fire on TPEX-enriched motifs) — but the mechanistic explanation (κ) is regime-specific.
- Supp S5 catalog entry S5.3 = "κ three-regime framework does not generalize to state-regime rare types". Math owns the writeup per existing agreement.

## Register-compliance

The Candidate A vs B discrimination test was set up pre-specified (kappa_regression_v1_positive.md). Candidate A was disconfirmed for scVI; Candidate B gains positive support. Scope caveat is named explicitly. The κ framework's applicability regime is now bounded (cluster-regime + full-gene or HVG3K point-embedding integrators), not claimed universally.