# Mathematics agent deliverables — master index

*agent-mozusggu, agent/mathematics, 2026-05-10. Head: `70bc6f6`.*

One-stop reference for the Mathematics agent's contributions to the CNS submission. Use for the author-contributions paragraph, reviewer cover-letter index, and to check that every theorem referenced in the manuscript has a home file.

## Theorems and lemmas

| Label | Name | File | Status |
|---|---|---|---|
| `thm:labelblind` / `thm:S1-blindspot` | Label-based blindness (Theorem 0) | `workspace/handoffs/abstract_thm0.tex`, `supplement_S1.tex` | Proved |
| `thm:minimax` / `thm:S1-minimax` | Minimax lower bound (i.i.d. + batch-het) | `supplement_S1.tex` §S1-minimax | Proved |
| `thm:S1-pt-stability` | P_t topological stability (Cohen-Steiner + DTM) | `supplement_S1.tex` §S1-pt-stability | Proved |
| `thm:S1-identifiability` | Cross-batch-witness identifiability | `supplement_S1.tex` §S1-identifiability | Proof sketch |
| `lem:exchangeability` / `lem:S1-exchangeability` | Exchangeability bound (Route C, bound form) | `supplement_S1.tex` §S1-exchangeability | Proved |
| `thm:S1-nulls` | Three-null calibration (N1/N2/N3) | `supplement_S1.tex` §S1-nulls | Proof sketches |
| `prop:rs-identifiability` / `thm:S1-real-closure` | RareShield + scVI ELBO identifiability | `supplement_S1.tex` §S1-real-closure | Proved |
| `thm:colm-identifiability` / `thm:S1-colm` | CoLM detectability on admissible F | `supplement_S1.tex` §S1-colm | Proved |
| `thm:S1-lecam-constants` | Le Cam two-point with explicit canonical-pop constants | `supplement_S1.tex` §S1-lecam-constants | Proved |
| `thm:S1-hierarchical` | Hierarchical REAL detectability amplification | `supplement_S1.tex` §S1-hierarchical | Proved |
| `rem:S1-extraction-annotation` | Compartment-classifier annotation bound | `supplement_S1.tex` | Documented |

## LaTeX drop-in fragments for manuscript

| File | Target location in BioSci manuscript |
|---|---|
| `workspace/handoffs/abstract_thm0.tex` | Abstract (one-sentence Theorem 0) |
| `workspace/handoffs/methods_drop_in.tex` | sections/07_methods.tex — Math formalization |
| `workspace/handoffs/sec02_results_lrs.tex` | sections/02_results_lrs.tex — Theorems 1-3 + exchangeability |
| `workspace/handoffs/sec06_results_rareshield.tex` | sections/06_results_rareshield.tex — RareShield theory |
| `workspace/handoffs/supplement_S1.tex` | supplement/S1_math.tex — full proofs |
| `workspace/handoffs/supp_table_S1_detectability.tex` | supplement/S1_math.tex — 23-row detectability table |
| `workspace/handoffs/supp_note_S0_8traps.tex` | supplement/S0_epistemics.tex — 8-trap mathematical mapping |
| `workspace/handoffs/hierarchical_real.tex` | sections/07 or supplement — Hierarchical REAL |
| `workspace/handoffs/theorem_zero_and_null_variance.tex` | supplement — alternate form + null-variance |
| `workspace/handoffs/label_map.tex` | preamble — label aliases reference |

## Supporting code

| File | Purpose | Signature |
|---|---|---|
| `workspace/code/compute_l_rare.py` | L_rare loss for RareShield training | `compute_l_rare(z_pre, z_post, motif_ids, cfg) -> (loss, diag)` |
| `workspace/code/hierarchical_real.py` | Telescoping-detection wrapper | `run_hierarchical_real(adata, integrator_fn, compartment_key, cfg) -> REALResult` |
| `workspace/handoffs/colm_loss_snippet.py` | Standalone CoLM loss (reference) | `colm_loss(z_pre, z_post, batch_id, k, quantile_rare) -> scalar` |

## Figures / tables

| File | Purpose |
|---|---|
| `workspace/handoffs/fig_ed1_envelope.pdf/png` | Fig ED1 (π, Δ) envelope at Kang-atlas scale |
| `workspace/handoffs/figED_pid_envelope.py` | Fig ED1 generator |
| `workspace/handoffs/table_s1_envelope.csv` | 23-pop table in CSV |
| `workspace/handoffs/supp_table_S1_detectability.tex` | Supp Table S1 (LaTeX) |
| `workspace/handoffs/ED_TPR_by_population.pdf/png` | ED bar chart, TPR per population colored by tier |
| `workspace/handoffs/ED_TPR_by_population.py` | ED bar chart generator |

## Reference notes (markdown, for cross-discipline reading)

| File | Audience | Purpose |
|---|---|---|
| `workspace/PLAN.md` | Math (self) | F1–F7 pipeline |
| `workspace/formal_model.md` | All | Measure-theoretic setup + Thm 2.1 (= Thm 0) |
| `workspace/detectability_framework.md` | DS, CB | D1–D4 label-free detectors with sample complexity |
| `workspace/lrs_theory_spine.md` | DS, BioSci, ML | 5-theorem theory spine in DS notation |
| `workspace/colm_admissibility.md` | ML | Admissible class F and CoLM detector |
| `workspace/handoffs/closed_form_tau.md` | Immunology | Closed-form τ + PH filtration levels |
| `workspace/handoffs/sub_compartment_pooling_bounds.md` | All | Hierarchical REAL math |
| `workspace/handoffs/predicted_detectability_registry_v1.md` | Immunology | Canonical predicted-detectability table |
| `workspace/handoffs/popid_crosswalk.md` | Immunology | Math ↔ registry_v1 ID map |

## Claims enabled in the manuscript

1. **Introduction / Abstract headline**: "Label-based evaluation is blind by construction to unannotated rare-subpopulation loss (Theorem 0)."
2. **Methods scalability claim**: "REAL detects populations with joint mass π and separation Δ whenever n·π²·Δ² ≥ 9 at α=β=0.05 on Kang 2024 (Λ=3, σ=1)."
3. **Methods identifiability claim**: "Cross-batch witnessing recovers a true rare motif with probability 1−δ when n_b·π_b² ≥ c_3·log(1/δ)/γ² per batch."
4. **Hierarchical REAL headline**: "Telescoping detection across compartment depth recovers populations up to 100× rarer than the whole-atlas floor permits."
5. **Exchangeability framing**: "The loss rate on held-out known rare motifs upper-bounds the expected loss rate on truly unannotated motifs of matched profile, up to η_match = L_F(η_π + η_Δ)."
6. **Honest disclosure**: "Of 22 rare target populations in the operating envelope, 20 are recovered (7 flat + 13 hierarchical); 2 are unresolved (DTP baseline, mast cells mode-4)."

## Author-contributions phrasing (draft for BioSci)

> "A.M. developed the mathematical framework: the label-based blindness theorem, minimax lower bound for label-free detection, topological and measure-theoretic stability results, identifiability proofs for REAL and RareShield, and the hierarchical-REAL detection-amplification theorem. A.M. also provided the LaTeX fragments and PyTorch modules that instantiate these claims in the manuscript and codebase."

---

*This index is living; updated with each commit. See `git log agent/mathematics` for history.*
