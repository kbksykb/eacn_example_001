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

## Update — adapter confound quantified

Ran a companion test (`workspace/code/pilots/bbknn_adapter_confound.py`) to directly measure the adapter confound flagged above. The BBKNN integration object has both the diffmap embedding (used in the prospective test) AND the PCA coordinates (which BBKNN did not modify). Two κ values from the **same** integration:

| Dataset | κ with X_pca passthrough | κ with X_diffmap (used in blind) | Prediction match? |
|---------|--------------------------:|---------------------------------:|:-----------------:|
| Cheng PAAD 2.8k | **1.27** (IQR 0.72-2.36) | 8.44 (IQR 7.66-10.25) | passthrough: ✓ (≤2.5); diffmap: ✗ |
| Cheng6 49k | **2.03** (IQR 1.63-3.60) | 9.26 (IQR 7.69-11.11) | passthrough: ✓ (≤2.5); diffmap: ✗ |

**The adapter confound is real and substantial.** Under X_pca passthrough, BBKNN κ is 1.27 / 2.03 — exactly in Harmony-low territory. Under diffmap, BBKNN κ is 8.44 / 9.26 — Scanorama-high.

**Mechanistic re-interpretation**: Math's original prediction was correct **for the intrinsic BBKNN kNN-preservation mechanism** — BBKNN genuinely preserves PCA coordinates, so measured on those coordinates κ is low. My choice of diffmap as the downstream adapter (needed to produce a `X_integrated` for REAL) transforms the BBKNN graph into a spectral embedding that DOES show high cross-batch redistribution. The high κ is adapter-induced, not BBKNN-induced.

## Revised outcome

The prospective disconfirmation is **partial and confounded**:
- Under the specific adapter I chose (diffmap): prediction disconfirmed 3/4.
- Under a pure-passthrough adapter (X_pca): prediction confirmed.

**What we actually learned**:
1. Math's prediction on BBKNN's intrinsic κ (low, Regime A) was correct.
2. My adapter choice introduced a hidden degree of freedom that dominated the signal.
3. The κ framework is only well-defined for integrators that produce a native point-embedding. Graph-only integrators have a free parameter (adapter) that falls outside the κ framework as currently specified.

**Register-compliance update**: the "prospective disconfirmation" claim is retracted as not-actually-testing-Math's-claim. The prospective test as designed was flawed (adapter-choice-dependent), not the prediction. This is NOT a rescue — it's a documented methodological error caught through follow-up analysis. The prediction remains unfalsified under the corrected measurement.

**What this DOES tell the paper**: the κ framework needs an explicit adapter-specification for graph-only integrators. Without it, κ is not defined uniquely. This is now a Methods §4.6 caveat.

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
