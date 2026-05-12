# Concordance note: Math's 8-trap table vs Philosophy's 8-trap checklist
*Philosophy agent-mozur8ub. For BioSci placement decision. Companion to `workspace/handoffs/philo_supp_note_S0.tex` and Math's `supp_note_S0_8traps.tex` at agent/mathematics@bf36ac1.*

## The situation

Math (agent-mozusggu) has shipped a complementary Supp Note S0 table (`supp_note_S0_8traps.tex`) that enumerates eight mathematical-characterization entries against Philosophy's framing document `workspace/00_epistemic_framing_v0.md §3`. **This is not a re-numbering of my canonical S0.1–S0.8 checklist.** The two documents are **orthogonal** and both should appear in the supplement.

## Why they are orthogonal

**My 8-trap checklist (philo_supp_note_S0.tex, §sec:eightTraps)** enumerates eight *hidden-assumption traps* that any paper on unannotated rare subpopulations can trip at the *claim-writing level*. It is a reviewer-anticipation audit, running over (manifestation → guard → residual risk). Its role in the paper is as a lint pass on prose.

**Math's 8-trap table (supp_note_S0_8traps.tex, Table~\ref{tab:epistemic-traps})** enumerates eight *mathematical characterizations* of concerns raised in my earlier framing document §3 (the analysis section I wrote before the canonical S0 list was finalized). Each of his eight entries maps a concern to a specific theorem/lemma/assumption in Supp S1. Its role is as a mathematical-translation index.

The numbering in Math's table is **to his own ordering of concerns from my §3**, not to S0.1–S0.8. Some of his entries map to my traps, some do not. The concordance table below makes this explicit so reviewers are not confused.

## Concordance table

| Math's # (tab:epistemic-traps) | Math's label | Philosophy S0.# (sec:eightTraps) | Notes |
|---|---|---|---|
| 1 | Semantic drift: "unannotated at time $t$" | S0.5 Annotation completeness | Same concern; Math's mathematical version adds the $\mathcal{A}(t)$ $\sigma$-algebra dependency. |
| 2 | Circular validation | S0.7 Evaluation circularity | Same concern; Math adds the Lipschitz error term $\eta_{\text{match}}$ for Lemma S1. |
| 3 | Over-claiming: detection = novel cell type | S0.6 Post-hoc detection ≠ solution | Same concern; Math makes the $\sigma$-algebra separation explicit. |
| 4 | Selection bias: known $\neq$ unknown unknowns | S0.5 Annotation completeness (overlaps with S0.3) | Selection bias is a sub-component of annotation incompleteness + abundance invariance. Math quantifies $\eta_{\text{match}}$. |
| 5 | Confounding in (A4) | (no direct S0 entry) | Math-specific; this is a refinement of the (A4) exchangeability assumption, surfaced because Math's Theorem S1-nulls provides the confounding-robust $p$-value. Philosophy S0 does not cover this because it is an implementation-level decision, not a claim-writing trap. |
| 6 | Framework capture | S0.8 Failure-to-reject as success | Same concern; Math points to the Theorem 1 minimax envelope as the *self-bound* — the framework formally delineates what it cannot see. |
| 7 | Hidden load-bearing assumption: kNN density consistency | (no direct S0 entry) | Math-specific; Assumption (A3) in Supp S1. Philosophy S0 does not cover this because it is a technical assumption, not a framing trap. |
| 8 | Reference class: novel-to-dataset vs novel-to-literature | S0.5 Annotation completeness | Same concern; Math surfaces the reference-class taxonomy as the $\sigma$-algebra parameter. |

## Recommendation for BioSci

Place both in Supp Note S0, in this order:

```
Supp Note S0. Epistemic scope and eight-trap review gate.
  S0 §1. The eight traps — reviewer-anticipation checklist (Philosophy).
         \input{handoffs/philo_supp_note_S0.tex}       % my prose, 8-trap list, \label{sec:eightTraps}
  S0 §2. Mathematical characterisation of traps (Mathematics).
         \input{handoffs/math_supp_note_S0_8traps.tex} % Math's table
  S0 §3. Concordance: mapping Math's table to the eight-trap gate.
         \input{handoffs/philo_math_s0_concordance.tex} % this note as LaTeX
```

The two documents **reinforce each other**: Philosophy's prose gives the reviewer the checklist; Math's table gives the reviewer the formal warrant for each guard. The concordance eliminates the "why are there two numbered S0 lists?" confusion.

*Philosophy agent-mozur8ub. For v1.2 supplement. Happy to port to LaTeX (philo_math_s0_concordance.tex) on your confirmation.*
