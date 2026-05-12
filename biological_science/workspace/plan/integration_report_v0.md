# Integration Report v0.2 — consolidating peer v0/v1 artifacts

Author: Biological Science Agent (`agent-mozuoaqx`)
Date: 2026-05-10
Status: internal coordination doc; input to manuscript v0.1

## What each peer has produced (as of branch fetch; updated v0.2)

| Discipline | Branch | Key artefact | Role |
|---|---|---|---|
| Data Science (agent-mozutvsy) | agent/data_science | `00_design_lrs_framework.md` (now REAL), `01_figure_plan.md` (6 mains + 10 ED), `02_lrs_benchmark_harness_spec.md`, `workspace/lrs_framework/`, `workspace/manuscript/figures/fig1/fig1.pdf` (delivered) | Detection framework lead; all figures |
| Mathematics (agent-mozusggu) | agent/mathematics | `formal_model.md` (Thm 2.1 label-blindness no-free-lunch), `detectability_framework.md` (4 label-free detectors + minimax sketch), `PLAN.md` | Methods Supp S1 + Thm 1 in main text |
| Immunology (agent-mozus3vv) | agent/immunology | `immunology_domain_brief.md`, `rare_immune_subsets_catalog.md`, `marker_panels.md` (oracle), `lrs_calibration_immunology.md`, `handoffs/05_results_kang_atlas_immuno_section.tex`, `handoffs/bib_additions_immunology.bib` | Intro immuno paragraph, Supp S4, Fig 3/4/6 interpretation |
| Tumor Biology (agent-mozurh1r) | agent/tumor_biology | `01_domain_brief_rare_subpop_oncology.md`, `02_rare_cancer_subpop_registry_v0.1.md`, `03_latex_fragments_for_biosci.tex`, `handoffs/tumor_bio_clinical.md` | Intro tumor paragraph, F4d/F4e cancer-type calls, Discussion clinical implications |
| Philosophy (agent-mozur8ub) | agent/philosophy | `00_epistemic_framing_v0.md` (A+B+C), `01_monitoring_log.md` | Intro paragraph 3, Discussion limitations, Supp Note S0 epistemic scope, 8-trap audit |
| Machine Learning (agent-mozurrmd) | agent/machine_learning | `literature/ml_synthesis.md`, `proposals/RareShield.md` (v2, vocabulary unified) | RareShield algorithm lead, Methods subsection, Intro/Discussion ML positioning |
| Computational Biology (agent-mozur9ik) | agent/computational_biology | `notes/experimental_plan_v1.md`, `manuscript/outline.md`, `code/rarescore/` prototype | Implementation + GPU runs, all benchmark outputs |

## Naming locked

| Layer | Name | Origin |
|---|---|---|
| Detection framework | **REAL** (Rare-subpopulation Erasure under Alignment, Label-free) | Consensus (my expansion); framework structure is DS's LRS spec |
| Per-seed score | **L** = HM(P_d, P_t, P_n, 1−B) | DS spec § 3 |
| Entity | **motif** (also: witnessed reproducible rare motif) | DS proposal |
| Motif set | **W** (witness set) | DS spec § 2 |
| Protection integrator | **RareShield** | ML agent |

Notation locked: R_b (rarity anchor), τ=2.0 (rarity threshold), σ_thresh=0.85 (cross-batch similarity), k=15, K=300, μ=2 (witnessing depth), β_w and λ_s (RareShield hyperparameters).

## Philosophy A+B+C as load-bearing frame

- **A** — population-level statistical realism about the unseen (Chao, Good-Turing) → informs uncertainty accounting.
- **B** — invariant/topological realism (what REAL does) → P_t, witness complex; Math formalises.
- **C** — exchangeability between held-out-known and truly-unknown rare subpopulations → calibrates L from positive controls and transfers as an upper bound; Math formalises.

Defensible strong claim (abstract's closing sentence + cover letter's central claim):

> We provide (i) a label-free, invariant-based detectability framework whose sensitivity is calibrated on held-out known rare subpopulations, and (ii) an integration strategy that demonstrably reduces the loss rate of such held-out rare subpopulations at pan-cancer-atlas scale, without degrading major-cell-type batch correction. We formalize the exchangeability assumption under which this loss rate transfers as an upper bound on loss of *unknown* rare subpopulations, and discuss the boundary conditions under which that assumption breaks.

Over-claiming past this is forbidden. Philosophy's 8-trap checklist is the pre-merge gate.

## Figure structure (6 mains + ≤10 ED per DS plan)

- F1 Conceptual (DS fig1.pdf landed v0); F2 REAL identifiability in silico; F3 Blinded pancreas/lung/retina; F4 Kang 4.9M; F5 RareShield ablation; F6 Discovery.

## Division of labor (unchanged from v0)

See first integration report; LaTeX fragment handoffs now exist from Tumor Biology and Immunology.

## Status summary

- **Manuscript skeleton**: v0.2 on agent/biological_science; sections 01–09 all have v0.1 prose.
- **Peer handoffs spliced**: Tumor-Bio Fig-4 + Discussion fragment; Immunology Kang subsection.
- **Open blockers**: Math LaTeX handoffs (Thm 1 + 5 theorems); ML Methods § RareShield; Philosophy Limitations subsection; Immunology Intro paragraph; Tumor-Bio Intro paragraph.
- **No CompBio results yet**: expected incrementally — Tier A (pancreas epsilon + HLCA ionocyte) first.

## Next actions

1. Commit v0.3 with REAL rename + peer handoff splicing. Push.
2. Submit result on t-mozv8bgw Milestone 1 (skeleton landed, ready to receive peer content).
3. Poll for Math LaTeX handoffs; when landed, splice into sections/02 + supplement.
4. Poll for ML RareShield Methods draft; splice into sections/09 § RareShield.
5. Poll for Philosophy Limitations draft; splice into sections/08.
6. Poll for Immunology Intro paragraph; splice into sections/01.
7. Poll for Comp-Bio first H5AD outputs; message DS to plot.
