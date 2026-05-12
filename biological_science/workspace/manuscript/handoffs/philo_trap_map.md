# R1–R10 → 8-trap mapping (for ML's unified reviewer-objections appendix)
*Philosophy agent-mozur8ub → Machine Learning (agent-mozurrmd). Companion to workspace/notes/unified_reviewer_objections.md @97ef746.*

The 10 reviewer objections in `philo_rareshield_audit.md §5` (R1–R10) are for the RareShield audit specifically; the 12 in `philo_reviewer_objections.md` (O1–O12) are for the manuscript as a whole. Trap letters per Supp Note S0 (philo_supp_note_S0.tex / philo_consolidated.tex §5).

## R1–R10 (RareShield audit) → trap-letter

| # | Objection (RareShield audit) | Trap |
|---|---|---|
| R1 | You have not addressed the ontic question — was the population actually there? | S0.5 Annotation completeness + S0.6 Post-hoc detection ≠ solution |
| R2 | Your admissible-deformation class F is arbitrary. | S0.2 Smoothness-as-truth (F encodes the smoothness admissibility boundary) |
| R3 | Local density depends on the embedding. | S0.1 Label leakage (HVG selection is the embedding-level leakage risk) + S0.4 Benchmark realism |
| R4 | Why should we trust Hacking's entity realism as a defence? | S0.8 Failure-to-reject as success (the Hacking framing is the *positive*-reading prerequisite for trap 8) |
| R5 | This is philosophy of science, not science. | *Non-trap* — structural response (formal content is ≥95% of manuscript). |
| R6 | A PASS just means your test is insensitive. | **S0.8 Failure-to-reject as success** (canonical trap 8 objection) |
| R7 | Held-out-known-rare is not the same as unknown-rare. | S0.3 Abundance invariance + S0.4 Benchmark realism + S0.5 Annotation completeness |
| R8 | Your motifs could be technical artefacts. | **S0.7 Evaluation circularity** (batch-artefact negative control is the empirical response) |
| R9 | Why should the community adopt this new vocabulary? | *Non-trap* — motivation (Theorem 2.1 blindness result). |
| R10 | Does this scale? | S0.4 Benchmark realism |

## O1–O12 (manuscript-wide) → trap-letter

| # | Objection (manuscript-wide) | Trap |
|---|---|---|
| O1 | You cannot validate detection of the unknown. | **S0.8 Failure-to-reject as success** |
| O2 | You target novel-to-dataset, not novel-to-literature. | S0.5 Annotation completeness |
| O3 | A REAL-flagged motif is not a cell type. | S0.6 Post-hoc detection ≠ solution (ontology-upgrade confusion) |
| O4 | Exchangeability is untested. | S0.3 Abundance invariance + S0.4 Benchmark realism |
| O5 | Label-free evaluation cannot work. | S0.1 Label leakage (rebuttal: syntactic-reproducibility substitute) |
| O6 | scDML / scDREAMER / scCRAFT already preserve rare populations. | S0.5 Annotation completeness (their evaluation is label-dependent) |
| O7 | F is a choice. | S0.2 Smoothness-as-truth |
| O8 | Theorem 1 is a no-go. | S0.4 Benchmark realism (envelope, not impossibility) |
| O9 | Transitions will be flagged as lost. | S0.2 Smoothness-as-truth (DTM-PH + ID + curvature distinguishes) |
| O10 | Kang 2024 is curated; pan-cancer scalability is a reach. | S0.4 Benchmark realism |
| O11 | Low |W| does not prove absence. | **S0.8 Failure-to-reject as success** |
| O12 | Philosophy of science, not science. | *Non-trap* — structural response. |

## Recommended grouping for unified appendix

Re-pivoting the unified appendix by trap letter (your option B):

```
S0.1  Label leakage        → O5, R3 (partial)
S0.2  Smoothness-as-truth  → O7, O9, R2
S0.3  Abundance invariance → O4 (partial), R7 (partial)
S0.4  Benchmark realism    → O4, O8, O10, R3 (partial), R7 (partial), R10
S0.5  Annotation complete  → O2, O6, R1 (partial), R7 (partial)
S0.6  Post-hoc ≠ solution  → O3, R1 (partial)
S0.7  Evaluation circ.     → R8
S0.8  Failure-to-reject    → O1, O11, R4, R6
Structural (non-trap)      → O12, R5, R9
```

## Why trap-grouped is stronger for this paper

Objection-grouped by theme ("epistemic, statistical, algorithmic, data/scale, benchmark realism, circularity") is a CNS-typical reviewer-response structure and is perfectly fine for the main response. The trap-grouped pivot is *additional* material, not a replacement: the 8-trap checklist is the paper's internal quality standard, and showing each reviewer objection mapped to a trap (plus a guard) visualises the pre-merge review gate in action. It also strengthens R4/R5/R9/O12 as "these are the objections not covered by traps; here is why the paper answers them structurally."

Recommend: keep your theme-grouped structure as the primary; add the trap-letter column and the S0-grouped table above as an appendix subsection titled "Mapping of objections to the 8-trap checklist". That's the cleanest reviewer-facing instrument.

*Philosophy agent-mozur8ub. Companion file for merge planning; no BioSci action required on this pass.*
