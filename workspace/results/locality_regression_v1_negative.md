# ℓ-locality compartment-amplification test — 9-point result (Math's falsifiable prediction)

**Date:** 2026-05-11
**Input:** agent/mathematics prediction on 9-point method×scale grid (LAMP3+ mregDC × Cheng cohorts × Harmony/Scanorama/scVI)
**Status:** **Negative result — ℓ does not collapse the method×scale grid onto a single scaling curve.** Honest Popperian outcome; reported per Philosophy register correction.

## Math's prediction (restated)

For each (method, dataset) point: define ℓ(method) := median over cells of |kNN(c, X_integrated) ∩ same-lineage(c)| / k, with k=30. Predict: `ℓ(method) · -log10(p_bh)` monotone in `n·π²·Δ²` across all points, and the 3 methods collapse onto one scaling curve when weighted by ℓ.

## Empirical result (9-point grid, k=30)

| Dataset | Method | n·π²·Δ² | p_bh | -log10(p_bh) | ℓ_coarse | ℓ_fine |
|---------|:------:|--------:|-----:|-------------:|---------:|-------:|
| PAAD 2.8k | Harmony | 0.21 | 1.000 | 0.00 | 1.000 | 0.967 |
| Cheng5 20k | Harmony | 2.93 | 0.0130 | 1.89 | 0.967 | 0.600 |
| Cheng6 49k | Harmony | 4.17 | 0.0041 | 2.39 | 0.967 | 0.667 |
| PAAD 2.8k | Scanorama | 0.21 | 0.005 | 2.30 | 1.000 | 1.000 |
| Cheng5 20k | Scanorama | 2.93 | 0.008 | 2.10 | 0.967 | 0.867 |
| Cheng6 49k | Scanorama | 4.17 | 0.0039 | 2.41 | 1.000 | 0.900 |
| PAAD 2.8k | scVI | 0.21 | 0.003 | 2.52 | 1.000 | 0.933 |
| Cheng5 20k | scVI | 2.93 | 0.004 | 2.40 | 0.967 | 0.800 |
| Cheng6 49k | scVI | 4.17 | 0.0030 | 2.52 | 0.933 | 0.733 |

Plot: `workspace/results/locality_grid_6point.png` (regenerated with 9 points)
Data: `workspace/results/locality_grid.csv`
Panel A: `workspace/results/panel_A_tight_v2.csv` (now 27 rows)

## Three distinct integrator behaviors

The 9-point grid reveals **three qualitatively different** detection regimes, NOT a single scale-monotonic curve:

1. **Harmony — scale-monotonic** (-log10(p_bh): 0 → 1.89 → 2.39). Preserves at small n, fires progressively at larger n. Matches Theorem 1 + §rem:S1-effective-floor-locality.

2. **Scanorama — mild scale-modulation** (-log10(p_bh): 2.30 → 2.10 → 2.41). Fires at all scales; weak non-monotone trend with scale. Integrator-intrinsic over-mixing regardless of n.

3. **scVI — scale-flat** (-log10(p_bh): 2.52 → 2.40 → 2.52). Fires uniformly at all three scales. Non-linear latent compression dissolves LAMP3+ motif independently of sample size.

## Disconfirmation of the ℓ-unification prediction

**Coarse ℓ** (DC/Mac/Mono/Mast/Neut): 0.93-1.00 across all 9 points — essentially invariant. Cannot discriminate between methods or predict -log10(p_bh) variation.

**Fine ℓ** (MajorCluster): differentiates methods (Harmony 0.60-0.97, Scanorama 0.87-1.00, scVI 0.73-0.93). But the ℓ·-log10(p_bh) points do NOT collapse into a monotonic ordering in n·π²·Δ². Counter-examples:
- Scanorama PAAD: n·π²·Δ²=0.21, ℓ_fine·-log10(p)=2.30
- Harmony Cheng5: n·π²·Δ²=2.93, ℓ_fine·-log10(p)=1.13
- scVI PAAD: n·π²·Δ²=0.21, ℓ_fine·-log10(p)=2.35

At a n·π²·Δ² **13× larger**, Harmony Cheng5 gives ℓ·signal HALF that of Scanorama PAAD. A monotone collapse would require opposite ordering.

## Additional caveat — scVI non-specificity

scVI on Cheng6 fires on 64.7% of non-rare motifs (non_rare_flag_bh=0.647), compared to Harmony (0.47) and Scanorama (0.50). scVI's aggressive firing is accompanied by reduced specificity — not pure detection gain.

## Implications

1. **Theorem 1 scaling** is empirically confirmed for Harmony specifically, not universally. Methods §4.6 now distinguishes "Regime A: scale-induced (Harmony)" from "Regime B: integrator-intrinsic (Scanorama+scVI)". Manuscript commit `a6f3114`.

2. **The ℓ-weighted unifying curve prediction** (Math's hypothesis) is **falsified** at 9 points. §rem:S1-effective-floor-locality **as an intra-method mechanism** (why Harmony fires below whole-atlas floor) survives. The cross-method ℓ unification claim does not.

3. **Supp S5 catalog update needed**: add "unified ℓ-weighted scaling curve" under "mechanisms considered and falsified" section, with this writeup as the supporting evidence.

4. **Fig 2 implications for Data Science**: consider ordering the 9 points by n·π²·Δ² but rendering method-level sub-scatter with different markers rather than expecting a single curve. The story is "method-dependent detection regimes", not "universal sigmoid detection scaling".

## Philosophy register compliance

This is a documented empirical disconfirmation of a mid-level falsifiable prediction. NOT reframed as method-ranking evidence. NOT reframed as partial support. Reported as negative result.

- The adjacent mechanism (compartment-locality per §rem:S1-effective-floor-locality) survives in a narrower scope: explains intra-method below-floor detection for Harmony specifically.
- The cross-method unification claim is retracted; Supp S5 catalog updated accordingly.
- scVI's high non-specificity (0.647 non-rare flag rate on Cheng6) is disclosed alongside its detection strength, NOT hidden behind the p_bh headline.
