# ℓ-locality compartment-amplification test — 6-point result (Math's falsifiable prediction)

**Date:** 2026-05-11
**Input:** agent/mathematics prediction on 9-point method×scale grid (awaiting scVI Cheng6 → reported at 6/9 to honor register discipline)
**Status:** **Negative result under both coarse and fine lineage definitions.** Honest Popperian outcome; reported per Philosophy register correction.

## Math's prediction (restated)

For each (method, dataset) point in the LAMP3+ mregDC × Cheng cohort scale-series:
- define ℓ(method) := median over cells of |kNN(c, X_integrated) ∩ same-lineage(c)| / k
- predict: `ℓ(method) · -log10(p_bh)` monotone in `n·π²·Δ²` across ALL points

Stated intuition: when compartment-locality differs across methods, weighting by ℓ should collapse the three methods onto a single scaling curve.

## Empirical result (6-point grid, k=30)

| Dataset | Method | n·π²·Δ² | -log10(p_bh) | ℓ_coarse | ℓ_fine | ℓ_coarse · -log10(p_bh) | ℓ_fine · -log10(p_bh) |
|---------|------:|--------:|-------------:|---------:|-------:|------------------------:|----------------------:|
| PAAD 2.8k | Harmony | 0.21 | 0.00 | 1.000 | 0.967 | 0.00 | 0.00 |
| Cheng5 20k | Harmony | 2.93 | 1.89 | 0.967 | 0.600 | 1.82 | 1.13 |
| Cheng6 49k | Harmony | 4.17 | 2.39 | 0.967 | 0.667 | 2.31 | 1.60 |
| PAAD 2.8k | Scanorama | 0.21 | 2.30 | 1.000 | 1.000 | 2.30 | 2.30 |
| Cheng5 20k | Scanorama | 2.93 | 2.10 | 0.967 | 0.867 | 2.03 | 1.82 |
| Cheng6 49k | Scanorama | 4.17 | 2.41 | 1.000 | 0.900 | 2.41 | 2.16 |

(lineage: DC / Mac / Mono / Mast / Neut at coarse; MajorCluster at fine)

Plot: `workspace/results/locality_grid_6point.png`

## Disconfirmation of the locality prediction

**Coarse ℓ (DC / Mac / Mono / Mast)**: ℓ≈0.97-1.00 for both methods at all three scales. No method differentiation; the `ℓ · -log10(p_bh)` points do NOT collapse into monotonic ordering in `n·π²·Δ²`. Scanorama PAAD (n·π²·Δ²=0.21, -log10(p)=2.30) sits **higher** than Harmony Cheng5 (n·π²·Δ²=2.93, -log10(p)=1.89) — a clear inversion, not a collapse.

**Fine ℓ (MajorCluster)**: ℓ differentiates the two methods (Harmony: 0.60-0.67 on Cheng5/6; Scanorama: 0.87-1.00). But the ordering still doesn't collapse — Scanorama PAAD at n·π²·Δ²=0.21 retains `ℓ_fine · -log10(p) = 2.30`, still above Harmony Cheng5's `1.13`.

**Mechanism-level reading**: Scanorama's strong signal on PAAD (2.8k) while Harmony preserves (p=1.0) is not explained by compartment-locality. Both methods preserve coarse lineage structure (DC/Mac/Mono stays discrete). What Scanorama does on 2.8k is **dissolve LAMP3+ cells into cDC2 at the intra-DC MajorCluster scale** (ℓ_fine=1.0 is misleadingly high — LAMP3+ cells get reassigned to *neighboring DC subtypes*, which are still "DC" in coarse and potentially still "M04_cDC2_CD1C"-proximal in fine).

## What this implies for the paper

1. **The compartment-amplification mechanism in its ℓ-weighted form does not survive empirical test at 6 points.** The "sub-compartment locality enables below-floor detection" narrative (§rem:S1-effective-floor-locality in the Math Supp) survives — `n_compartment ≈ n_myeloid ≫ n·π²·Δ²` IS why Cheng5/Cheng6 fire below whole-atlas floor. But the cross-method prediction via ℓ is falsified.

2. **The Method×Scale story needs to be stated carefully.** Current Methods §4.6 says "Theorem 1 scaling predicts real-data detectability" — still empirically supported for Harmony (monotone 0 → 1.89 → 2.39). But claiming Theorem 1 predicts Scanorama's behavior requires additional mechanism — Scanorama fires at PAAD (n·π²·Δ²=0.21, well below floor) where Theorem 1 predicts preservation. Scanorama's detection-at-2.8k is **integrator-aggressiveness**, not scale-induced compartment-amplification.

3. **Honest positioning for BioSci Methods**: Theorem 1 predicts Harmony scale-monotonicity. Scanorama's cross-scale firing pattern is method-intrinsic (always-aggressive), not scale-predicted. Two orthogonal mechanisms, both documented.

## Action items

- [ ] Update `workspace/manuscript/section_04_methods_draft.md` §4.6 to distinguish scale-induced (Harmony) from method-intrinsic (Scanorama) detection regimes.
- [ ] Ping Math with the 6-point result so the formal ℓ-statement gets revised in the Supp before Supp freeze. Math owns the mathematical treatment; I own the empirical disconfirmation.
- [ ] Once scVI Cheng6 lands, extend to 9 points. Prediction: scVI's ℓ sits between Harmony and Scanorama per Math's original estimate (ℓ_scVI≈0.85); 3rd-method signal may or may not fit the unifying curve.
- [ ] Consider a different unifying variable. Candidates: fraction-of-rare-cells-retained-within-same-MajorCluster (per-cell preservation of held-out LAMP3+ cells specifically), OR the OC-candidate count (how many non-rare motifs also flag).

## Philosophy register compliance

This is a documented empirical disconfirmation of a falsifiable mid-level prediction. NOT reframed as method-ranking evidence; NOT reframed as partial support. Called a negative result.
