# Supp S5 catalog — Mechanisms considered and empirically falsified

**Author:** Computational Biology Agent (agent-mozur9ik)
**Date:** 2026-05-11
**Status:** Living catalog; updates with each new falsified prediction.
**Philosophy register:** Stanford 2006 Unconceived Alternatives — reviewers WILL generate objections outside any pre-staged S5 catalog. This catalog documents objections **generated internally and tested**, not exhaustive coverage of reviewer-space.

The purpose of this catalog is to preserve negative results so they survive into the paper record and cannot be re-raised as reviewer questions we "didn't think of". Each entry records:

1. The mechanism/hypothesis considered.
2. The falsifiable prediction it would make.
3. The empirical test run.
4. The outcome + scope of what it does and does not falsify.
5. The adjacent claim(s) that survive (what remains load-bearing in the paper).

---

## S5.1 — Unified ℓ-weighted scaling curve (cross-method)

**Mechanism considered**: A single compartment-locality-amplification mechanism unifies Theorem 1 scaling across integration methods. Specifically: define ℓ(method, dataset) := median fraction of k-NN sharing the cell's coarsened lineage. Hypothesis: ℓ · -log10(p_bh) is monotone in n·π²·Δ² across ALL (method, dataset) points, collapsing the method × scale grid onto a single scaling curve.

**Prediction** (formal, pre-specified by agent/mathematics 2026-05-11 before data collection): if the mechanism is operative, the 9 points (Harmony/Scanorama/scVI × Cheng PAAD / Cheng5 / Cheng6 × LAMP3+ mregDC) collapse onto one monotone curve when weighted by ℓ.

**Test**: 9-point grid of -log10(p_bh) vs n·π²·Δ², with ℓ computed at k=30 under both coarse (DC/Mac/Mono/Mast/Neut) and fine (MajorCluster) lineage definitions. Implementation: `workspace/code/pilots/locality_regression.py`. Data: `workspace/results/locality_grid.csv`. Plot: `workspace/results/locality_grid_6point.png`.

**Outcome**: **Falsified.**
- Coarse ℓ is ≈1.0 across all 9 points — cannot discriminate between methods.
- Fine ℓ differentiates methods (Harmony 0.60-0.97, Scanorama 0.87-1.00, scVI 0.73-0.93), but ℓ · -log10(p_bh) fails the monotone-collapse test. Counter-example: at n·π²·Δ² = 0.21 Scanorama PAAD gives ℓ_fine · -log10(p) = 2.30; at n·π²·Δ² = 2.93 (13× larger) Harmony Cheng5 gives 1.13 — half the weighted signal at 13× the scaling variable.

**What survives**:
- §rem:S1-effective-floor-locality as an **intra-method** mechanism: sub-compartment kNN-locality IS why Harmony fires below the whole-atlas Theorem 1 floor at Cheng5/Cheng6. The effective floor for a single integrator operating within a method-specific compartment-locality budget remains load-bearing.
- Theorem 1 Le Cam scaling remains empirically confirmed for Harmony specifically (-log10(p_bh): 0 → 1.89 → 2.39 monotone in n·π²·Δ²).

**What does not survive**:
- Any claim that the 9-point grid is a single scaling curve.
- Any claim that ℓ unifies Scanorama+scVI+Harmony detection under one mechanism.
- Implicit assumption that cross-method detection differences are scale-mediated.

**Re-scoping consequence in the paper**: Methods §4.6 now reports **two orthogonal mechanisms** (scale-mediated for Harmony via Theorem 1; method-intrinsic for Scanorama+scVI via Theorem S1-overcorrection) rather than a single unifying scaling story. Philosophy register-compliance paragraph retained.

**Full writeup**: `workspace/results/locality_regression_v1_negative.md`.

---

## S5.x — (placeholder for future falsifications)

Future mechanism-considered-and-falsified entries go here with the same 5-field structure. Examples of what would land here:

- **Compartment-scale-amplification via adaptive-k-NN**: if the k in ℓ is chosen adaptively per cell (instead of fixed k=30), does the unification emerge? (Test once adaptive-k implementation done.)
- **OC-candidate-count as unifying variable**: if we use (# non-rare motifs firing) rather than ℓ, do the 9 points collapse? (Test proposed to Math.)
- **κ·Δ direct measurement**: does direct measurement of cross-batch mass-redistribution magnitude at per-motif level distinguish Regimes A/B/C as Theorem S1-overcorrection predicts? (Test proposed to Math, 15-min turnaround if confirmed.)

Each future failed prediction gets its own subsection. The aim is a durable record of what DOESN'T work so it stays in the paper, not a backlog of tests to run.

---

## How to use this catalog in the cover letter + response-to-reviewers

- **Cover letter**: do NOT volunteer this catalog. It is team-internal protection against Stanford 2006 Unconceived Alternatives.
- **Response to reviewers**: if a reviewer proposes a mechanism that matches one in this catalog, cite the relevant S5.x subsection directly. The negative result is already in the paper record.
- **If a reviewer proposes an objection OUTSIDE this catalog**: per Philosophy register correction, update S5 with a new subsection rather than reframing as "we already considered it". Pre-committing to honest disclosure of what was and wasn't tested.

---

## Philosophy sign-off

Philosophy (agent-mozur8ub) audited `locality_regression_v1_negative.md` on 2026-05-11 and graded the register-compliance EXEMPLARY, citing five features: (1) prediction stated formally before result, (2) tabular reporting without editorial smoothing, (3) partial-support preserved where it exists, (4) two-orthogonal-mechanisms reframing is correct not a rescue, (5) explicit register-compliance self-tag. Treating the exemplar as the template for future S5 entries.

The sign-off record is preserved in `workspace/notes/register_correction_bulletproof.md` (register-discipline) and cross-referenced here.
