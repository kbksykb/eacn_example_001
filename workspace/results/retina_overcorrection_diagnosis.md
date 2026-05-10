# Retina overcorrection motif diagnosis — RBC sub-fragment redistribution, not rare-loss

**Commit:** agent/computational_biology@[latest]
**Date:** 2026-05-11

## What was flagged

On Shekhar 2016 retina (BC8_9, BC4, BC2, BC1A holdouts × Scanorama + scVI = 8 integrations), REAL's OT channel consistently flagged 4-5 motifs as "loss" with p<0.05 and loss_probability > 0.65. These were ALL non-rare-truth motifs by the holdout's ground-truth labeling.

## What the motifs actually are

Mapping Leiden over-cluster IDs to Shekhar annotations:

| Motif | n cells | Top Shekhar label | % pure |
|------:|--------:|-------------------|-------:|
|     5 |     908 | RBC               | 99.8% |
|     7 |   2,197 | RBC               | 99.6% |
|     8 |   1,007 | RBC               | 100%  |
|     9 |   2,171 | RBC               | 98.8% |
|    11 |     638 | RBC               | 97.0% |
|    16 |     130 | RBC               | 99.2% |
|    21 |      43 | RBC               | 100%  |

**RBC = Rod Bipolar Cells**, the most abundant Shekhar type (8,175 cells = 41% of the dataset). Leiden over-cluster at resolution=2 fragmented the single biological RBC population into **5-7 sub-motifs** due to internal transcriptional heterogeneity visible in the pre-integration 50-dim PCA.

## What Scanorama/scVI do

When these integrators align the batches, they collapse the 5-7 RBC sub-fragments toward a single tight post-integration cluster. This redistributes mass: cells previously in 7 distinct density modes end up in 1. REAL's CoLM channel correctly detects this as mass redistribution violating the admissible batch-effect field class F.

Harmony, being more conservative, does not collapse the sub-fragments as aggressively and therefore does not trigger the flag.

## Why this is NOT what we originally framed

Earlier messages called these motifs "overcorrected abundants" with the implication that REAL was detecting a separate failure mode (generic overcorrection of abundant structure). That framing is **partially right but imprecise**:

- **Correct**: REAL's CoLM fires on any mass redistribution violating F. These are real redistribution events.
- **Imprecise**: calling them "overcorrection of abundant" suggests the biology is damaged. What's happening is that Leiden's pre-integration over-clustering artefactually fragmented RBC, and Scanorama/scVI collapse the fragments back to biological coherence. From a pure biology standpoint, the RBC re-unification is arguably CORRECT integration behaviour, not overcorrection.

## Reframed interpretation for the paper

REAL detects **three distinguishable regimes** on single-cell batch integration:

1. **Rare-subpop loss** (synth factorial α=0.4 regime): detector fires on a biologically distinct rare population that got absorbed into an abundant one. The Fig 2b synthetic story.

2. **Abundant-sub-structure redistribution** (retina RBC regime): detector fires on mass movement between fragments of an abundant population. This is NOT necessarily "overcorrection" — it may reflect the integrator doing its job of unifying artefactually-split clusters. Whether this is "good" or "bad" depends on downstream biology.

3. **Pure batch-artefact** (E_neg regime): detector does NOT fire (OT p=1.000) even when a batch-specific tight density mode looks superficially like a rare subpop. Fig 3d specificity.

**Manuscript implication**: Fig 2c / Fig 3d framing must distinguish regime 2 from regime 1. The retina experiment shows REAL is correctly sensitive to regime 2 (mass movement), but regime 2 is not "rare-subpop loss" — it's pre-integration over-cluster fragment redistribution. The Fig 2c title "retina specificity = TRUE NEGATIVE on BC8/9 rare-truth" captures this properly.

## Downstream actions

1. Update workspace/manuscript/section_03_draft.md to distinguish regime 2 from regime 1.
2. DS's Fig 2c already positioned as specificity/true-negative panel per latest ingest (agent/data_science@ef6b678). That is correct given this diagnosis.
3. For RareShield A/B on retina: the "Δlp -0.236" result is still meaningful because RareShield DID reduce the CoLM signal on these motifs. But the CLAIM should be "RareShield reduces abundant-sub-fragment-redistribution signal without harming major-type ARI" — NOT "RareShield rescues rare-subpop loss". These are different biology claims with different reviewer-relevance.
4. Philosophy's trap-7/8 analysis is STRONGER now: REAL correctly flags regime 2 (real structural redistribution) and does NOT flag regime 3 (pure artefact). That's a clean 2×2 truth table.

## Updated evidence pack

**Rare-subpop loss detection (the headline claim)**: synthetic factorial at α=0.4, Scanorama flags rare-truth at p=0.0099, abundants not flagged at p≤0.32. **This remains the only controlled demonstration of Regime 1.**

**Abundant-sub-fragment redistribution detection (secondary claim)**: retina RBC sub-fragments under Scanorama/scVI. Novel observation, but semantically different from rare-subpop loss.

**Specificity / pure artefact**: E_neg synth, both Harmony and Scanorama give OT p=1.000.

**Protection**: RareShield A/B drops loss-probability on flagged motifs while keeping major-type ARI within gate. Works on both Regime 1 (synth) and Regime 2 (retina).
