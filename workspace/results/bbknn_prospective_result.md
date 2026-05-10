# BBKNN prospective prediction — result: DISCONFIRMED

**Date:** 2026-05-11
**Pre-registration:** `workspace/results/bbknn_prospective_prediction.md` (Math's formal prediction, committed before data observation at commit da01cac)
**Status:** BLIND RUN COMPLETE. Results inspected. **3 of 4 pre-registered criteria failed.**

## Result numbers

| Dataset | κ_median | IQR | rare_min_p_raw | rare_min_p_bh | non_rare_flag_bh | #rare_motifs |
|---------|---------:|----:|---------------:|--------------:|-----------------:|------------:|
| Cheng PAAD 2.8k | **9.228** | [8.472, 11.424] | 0.0020 | **0.0105** | 0.200 | 1 |
| Cheng6 49k | **10.335** | [8.680, 13.032] | 0.0020 | **0.0044** | 0.441 | 1 |

## Prediction evaluation

| Criterion | Prediction | Observation | Pass? |
|-----------|------------|-------------|:-----:|
| κ_med Cheng PAAD ≤ 2.5 | ≤ 2.5 | **9.228** | ✗ |
| κ_med Cheng6 ≤ 2.5 | ≤ 2.5 | **10.335** | ✗ |
| Cheng PAAD p_bh ≥ 0.5 (preserve) | ≥ 0.5 | **0.0105** | ✗ |
| Cheng6 p_bh < 0.01 (fire) | < 0.01 | 0.0044 | ✓ |

**Outcome: DISCONFIRMATION** (3/4 criteria failed).

## What was disconfirmed

Math's hypothesis that "BBKNN preserves kNN structure → low-κ → Regime A" was disconfirmed at both scales. BBKNN's κ (9.2-10.3) lies in Scanorama-high territory (κ 8.3-9.1 on Cheng PAAD/5/6), not Harmony-low territory (κ 2.09-2.35). BBKNN fires at Cheng PAAD 2.8k (p_bh=0.01), matching Scanorama's method-intrinsic-at-any-scale behavior, not Harmony's preserve-at-small-n.

## Important confounding disclosure — adapter choice

The BBKNN integration adapter used here emits `X_integrated = X_diffmap(n_comps=32)` computed on the BBKNN kNN graph (scanpy.tl.diffmap on the bbknn-modified connectivities). BBKNN does NOT produce a native point-embedding — it produces a batch-balanced kNN graph. The adapter choice matters:

- **diffmap embedding** (our choice): picks up spectral structure of the BBKNN graph, which can amplify cross-batch signal if the graph flattens local batch clusters.
- **Alternative: use X_pca directly** (BBKNN doesn't change PCA coords, only neighbors): would give κ≈0 trivially because pre and post X_pca are identical.
- **Alternative: UMAP embedding** on BBKNN graph: a different spectral reduction.

This is a methodological limitation of the test. Math's prediction framed BBKNN as a "kNN-preserving integrator", which implicitly assumes κ should be measured on the *embedding* that downstream REAL scoring actually uses. Our diffmap adapter is the reasonable choice for an honest REAL pipeline (diffmap is what users would actually use for visualization/clustering after BBKNN), but other adapters could yield different κ.

**What this means for the paper**: the κ three-regime framework assumes a point-embedding `X_integrated` exists. For graph-only methods like BBKNN, the adapter choice is a hidden degree of freedom. This is a CAVEAT to scope explicitly in Methods §4.6, not a rescue attempt for the disconfirmed prediction.

## Mechanistic re-interpretation (post-disconfirmation)

Given BBKNN + diffmap κ ≈ 9-10 (Scanorama-like) and BBKNN-on-Cheng-PAAD p_bh = 0.01 (fire):
- BBKNN + diffmap appears to land in **Regime B (method-intrinsic), not Regime A** as predicted.
- The batch-balanced kNN graph, when converted to a spectral embedding, mixes across batches more aggressively than Harmony's linear correction — consistent with the "BBKNN is aggressive mixing" intuition some practitioners hold, and AGAINST Math's "kNN preservation" intuition.
- The three-regime taxonomy survives as descriptive categories, but **the method → regime mapping is not predictable from design intuition** — it must be measured empirically.

## What this does and does not falsify

**Falsifies**: Math's specific prediction that "BBKNN-preserves-kNN" implies low-κ. The mechanistic interpretation of κ as "violation of kNN structure" is too narrow.

**Does NOT falsify**:
- The three-regime taxonomy itself (low / medium / high κ stratifies methods).
- The ℓ-unification negative result (still stands, separate experiment).
- The intra-method compartment-locality mechanism (still stands for Harmony specifically).

## Register-compliance notes

- Prediction was registered BEFORE run (commit da01cac).
- Runner was blind until this evaluation file was written (commit after).
- Disconfirmation conditions were pre-specified in the registration.
- Result is reported at face value; no re-fitting of criteria post-hoc.
- Adapter-choice confound is disclosed in the body, not buried.

**Three falsifiability tests now in the paper record**:
1. ℓ-unification retrospective, falsified (Supp S5.1).
2. BBKNN-κ prospective, falsified (Supp S5.2 candidate).
3. Sade-Feldman κ three-regime, partial disconfirmation (pending HVG3K control).

Per Philosophy's observation: the paper has a documented habit of prediction-test-disconfirmation, which makes the un-disconfirmed claims (Theorem 1 Harmony scale-series, κ three-regime stratification on Cheng+Zheng cluster-regime, ΔARI gate on RareShield A/B) more credible.

## Next steps

- Supp S5.2 entry documenting this disconfirmation (Math owns Supp S5 per their message).
- Consider re-running BBKNN with X_pca-unchanged adapter to see if κ changes regime. If it does, the adapter-choice confound is real; if κ stays high, the BBKNN-mixes-aggressively interpretation wins.
- Methods §4.6 addendum: κ framework applies to point-embedding integrators; graph-only integrators (BBKNN, some kNN-based methods) require an adapter choice that is itself a κ-affecting variable.
