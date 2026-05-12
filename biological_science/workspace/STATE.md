# Biological Science Agent — State Snapshot

**Agent ID**: `agent-mozuoaqx`
**Branch**: `agent/biological_science`
**Latest commit**: (see `git log --oneline -1`)
**Manuscript state**: 29-page main + 10-page supplement + 3-page cover letter, PDF-compiled via pdflatex (TeX Live 2022).

## Compiled artefacts

- `workspace/manuscript/main.pdf` — main manuscript (29 pages)
- `workspace/supplement/supplement.pdf` — supplementary information (10 pages)
- `workspace/cover_letter/cover_letter.pdf` — cover letter to CNS editors (3 pages)

## Source files

- `workspace/manuscript/main.tex` — top-level main manuscript
- `workspace/manuscript/sections/01..09_*.tex` — nine content sections (intro → methods)
- `workspace/supplement/supplement.tex` — supplement with S0–S9 notes
- `workspace/cover_letter/cover_letter.tex` — cover letter
- `workspace/manuscript/figures/{fig1..fig5,ed1}/*.pdf` — vendored figure PDFs from agent/data_science
- `workspace/manuscript/handoffs/*.tex` — 47 peer LaTeX fragments, spliced by \input

## Peer contributions spliced

| Discipline | Role | Artefacts |
|---|---|---|
| Data Science (agent-mozutvsy) | Framework spec, all figures, Methods §5 | LRS framework, Figs 1-5 + ED1, figure_captions.md, labels.tex |
| Mathematics (agent-mozusggu) | Theorems + Supp S1 | Thm 0 (label-blind), Thm 1-3 (minimax / stability / identifiability), Lem S1 (exchangeability bound), Thm S1-hierarchical, 11 lemmas total |
| Machine Learning (agent-mozurrmd) | RareShield algorithm | Unified L_RS loss form, Methods §RareShield, lit positioning, unified reviewer-objections pack |
| Computational Biology (agent-mozur9ik) | Implementation + experiments | Pancreas-synth pilot, factorial sweep v1 (Scanorama p=0.0099, Harmony no-flag), RareScore → REAL channels code |
| Immunology (agent-mozus3vv) | Lineage/state, Supp S2-S3 | Intro paragraph, Kang stress-test rationale, Discussion implications, limitations, oracle-score Methods, below-floor close |
| Tumor Biology (agent-mozurh1r) | Clinical/TME, Discussion | Intro paragraph, Fig 4e cancer-type table, Fig 5 melanoma clinical framing, clinical-implications Discussion |
| Philosophy (agent-mozur8ub) | Epistemic framing | Measured framing, Meno Intro P3, onco-intro para, Limitations, Below-floor head, Supp S0 8-trap checklist, cover-letter philosophy paragraph, lint:epistemic tool |

## Pipeline health

- `lint:epistemic` → PASSES clean on manuscript/sections/ (no ungarded red-flag phrases).
- `pdflatex` → compiles cleanly on TeX Live 2022 with cleveref + longtable + enumitem + amsthm.

## Outstanding blockers (paper completion)

1. CompBio real-data H5AD (pancreas and/or retina). Network stack is unreliable from this box; human-stage may be needed.
2. Kang 2024 integration run → Table 1 motif IDs → Immuno + Tumor-Bio 24h SLA annotation → Fig 4e stacked bar finalized.
3. RareShield training on synth → Fig 3 ablation columns → Fig 3 full.
4. Merge all bib fragments (handoffs/immuno_bib.bib + Math bib + Philo bib) into single references.bib.

## Versions shipped

| Version | Commit | Highlights |
|---|---|---|
| v0.1 | 6caa67c | Scaffold CNS manuscript skeleton |
| v0.2 | 81023e3 | Unify LRS/RareShield; absorb peer seeds; 6-fig plan |
| v0.3 | af478de | REAL/RareShield rename; splice Tumor-Bio + Immuno handoffs |
| v0.4 | 62feb29 | Splice Math Theorems 0-5; Immuno Intro |
| v0.5 | d5aa1fa | Math v2 + ML RareShield + all peer v0.1 handoffs |
| v0.6 | 08b3c0c | ML unified RareShield; CompBio first pilot |
| v0.7–v0.8 | 4dff4e5 | All Philosophy handoffs; Supp S0; Thm 0 in abstract |
| v0.11 | 541e9f4 | Vendor DS figure PDFs; embed \includegraphics |
| v1.0-rc1 | 9b7d37e | **First PDF compile**: main 22pp + supp 8pp + cover 2pp |
| v1.1 | fbef79e | Philo §02+§08 power-floor stratification; Math hierarchical REAL; 27pp |
| v1.2 | 4117111 | Philo below_floor_head; CompBio factorial sweep; phantom \Cref; 28pp |
| v1.3 | e9f9757 | Philo onco-intro + Immuno below-floor-close + Math Supp Table S1; 28pp / 10pp |
| v1.4 | 3f1550e | Supp S9 reviewer objections; lint:epistemic passes |
| v1.5 | cd8f159 | Tumor Fig 5 melanoma clinical + Immuno cover-letter paragraphs; 29pp / 10pp / 3pp |
| v2.0 | a643313 | FIRST REAL-DATA HEADLINE: pancreas LossRate@10 0.88/0.98/1.00 + retina 0.17/0.43/0.44 |
| v2.1 | d74f2cf | Math Thm S1-overcorrection in §08; Philo Meno-template bookend |
| v2.2 | 7fd232d | Math "envelope binds in practice" paragraph; Immuno specificity reframing |
| v2.3 | 3f60ef5 | Philo v1.8 broadening: REAL as overcorrection-detection, rare-subpop as worst-case |
| v2.4 | e443a79 | Merge peer bibs → references.bib; plainnat style; Philo low-priority audits applied. **Main 33pp with rendered bibliography** |
| v2.5 | 9a16344 | CompBio RareShield A/B synth (AUPRC 0.556→0.867, +0.31, ΔARI=0); retina 9×3 abundance all at-floor per Thm 1 |
| v2.6 | 3855e71 | Philo philosophy-of-science bib backfill (20 entries); reviewer-preempt v2 (28 objections); §08 consolidated |
| v2.7 | 12b40f9 | Path-1 root-cause framing absorbed — Triple-Absence loop (§01); convergence-vs-dispersion (§02); circular-reasoning caveat (§08) |
| v2.8 | 20343cb | **Retina RareShield A/B on all 3 methods** — Scanorama Δlp=-0.236 on overcorrected abundants, ΔARI=-0.014, BC8_9 unchanged; Harmony control +0.004 ARI, scVI partial. **First real-data protection demonstration.** Main 34pp. |
