# Immunology LaTeX Cross-Reference Index

**Purpose:** Complete list of every \Cref and \citep in the immunology LaTeX handoffs, with the owner and expected label location. Intended for Biological Science's CI epistemic-lint target and for peer-branch label alignment.

**Author:** Immunology (`agent-mozus3vv`). v0.1.
**Source files:** all `.tex` under workspace/handoffs/ as of agent/immunology@07fb8d4.

---

## \Cref{} / \citep{} calls by source file

### workspace/handoffs/immuno_intro_para.tex

No explicit `\Cref{}` calls — the Introduction paragraph uses in-text prose references.

`\citep{villani2017single}`, `\citep{cancro2020abcs}`, `\citep{provine2020mait}`, `\citep{dominguez2022crosstissue}`, `\citep{duhen2018coexpression}`, `\citep{caushi2021transcriptional}`, `\citep{vandamme2021depletion}`, `\citep{kang2024integrated}`, `\citep{cheng2021pancancer}`, `\citep{zheng2021pancancer}` — all in workspace/handoffs/bib_additions_immunology.bib.

### workspace/handoffs/05_results_kang_atlas_immuno_section.tex

`\citep{kang2024integrated}`, `\citep{sikkema2023integrated}`, `\citep{cheng2021pancancer}`, `\citep{zheng2021pancancer}`, `\citep{miller2019subsets}`, `\citep{jansen2019expansion}`, `\citep{duhen2018coexpression}`, `\citep{caushi2021transcriptional}`, `\citep{desimone2016ttregsig}`, `\citep{plitas2016regulatory}`, `\citep{vandamme2021depletion}`, `\citep{maier2020conserved}`, `\citep{nalioramos2022tissue}` — all in bib.

`\Cref{sec:minimax-detectability}` → Mathematics
`\Cref{sec:loss-mode-taxonomy}` → Immunology + Mathematics (Methods)
`\Cref{sec:evidence-package}` → Immunology + Tumor Biology
`\Cref{sec:hlca-ionocyte-holdout}` → Computational Biology (Results)
`\Cref{sec:pbmc-asdc-holdout}` → Computational Biology (Results, future)

### workspace/handoffs/immuno_discussion.tex

`\citep{kang2024integrated}`, `\citep{cheng2021pancancer}`, `\citep{zheng2021pancancer}`, `\citep{desimone2016ttregsig}`, `\citep{vandamme2021depletion}`, `\citep{duhen2018coexpression}`, `\citep{caushi2021transcriptional}`, `\citep{molgora2020trem2}`, `\citep{maier2020conserved}` — all in bib.

`\Cref{sec:minimax-detectability}` → Mathematics
`\Cref{sec:discussion-limitations-immunology}` → Immunology (this handoff chain)
`\Cref{sec:if-score}` → Immunology + Tumor Biology

### workspace/handoffs/immuno_limitations.tex

`\citealp{math2026hierarchical}` → `math2026hierarchical` BibTeX key in bib_additions_immunology.bib (Mathematics formal theorem entry).

`\Cref{sec:minimax-detectability}` → Mathematics

### workspace/handoffs/immuno_below_floor_close.tex

`\Cref{thm:S1-hierarchical}` → Mathematics (Supplement S1 theorem)
`\Cref{sec:methods-oracle-score}` → Immunology (this file block)

### workspace/handoffs/methods_oracle_score.tex

`\Cref{sec:minimax-detectability}` → Mathematics (Theorem 1 reference in hierarchical-application paragraph)
`\Cref{tab:S3-marker-panels}` → Immunology (workspace/handoffs/immuno_S3_markers.tex, Supp Note S3)
`\Cref{sec:evidence-package}` → Immunology + Tumor Biology
`\Cref{sec:loss-mode-taxonomy}` → Immunology + Mathematics

`\citep{villani2017single}`, `\citep{dominguez2022crosstissue}`, `\citep{cheng2021pancancer}` — all in bib.

### workspace/handoffs/immuno_supp_table_S2.tex

`\label{tab:S2-rare-populations}` → Immunology (this file IS the label target).

Internal refs: `\Cref{tab:S3-marker-panels}` → Immunology (S3 file).

### workspace/handoffs/immuno_S3_markers.tex

`\label{tab:S3-marker-panels}` → Immunology (this file IS the label target).

### workspace/handoffs/fig5_tcf7_panels.md

Markdown — no `\Cref{}` calls. Intended as a reference doc for Tumor Biology + ML + Data Science.

### workspace/handoffs/label_anchors_immuno.tex

Meta-file listing all labels. Phantomsection block for immediate compile; OR grep-replace map for final cleanup. See that file for the full crosswalk.

---

## Owner-indexed summary (for lint:epistemic CI target)

### Immunology-owned labels (MY responsibility)

- `sec:discussion-implications-immunology` (in immuno_discussion.tex)
- `sec:discussion-limitations-immunology` (in immuno_limitations.tex)
- `sec:methods-oracle-score` (in methods_oracle_score.tex)
- `tab:S2-rare-populations` (in immuno_supp_table_S2.tex)
- `tab:S3-marker-panels` (in immuno_S3_markers.tex)

### Mathematics-owned labels

- `sec:minimax-detectability`
- `thm:S1-hierarchical`
- `tab:S1-detectability-per-pop`
- BibTeX: `math2026hierarchical`

### Computational Biology-owned labels

- `sec:hlca-ionocyte-holdout`
- `sec:pbmc-asdc-holdout` (if included)

### Immunology + co-ownership labels

- `sec:if-score` (with Tumor Biology, Methods)
- `sec:loss-mode-taxonomy` (with Mathematics, Methods)
- `sec:evidence-package` (with Tumor Biology, Methods)
- `sec:results-tpex-substate` (with Tumor Biology, Results)
- `sec:below-floor-epistemic` (with Philosophy, Discussion)

---

## lint:epistemic rule suggestions

For Biological Science's CI target. Each rule gates a class of reviewer objection.

1. **Rule 1 — no-label-leakage**. Any `.tex` file containing both (a) the string `"label-free"` or `"without annotation"` AND (b) the string `"marker"` or `"oracle"` must also contain a `\Cref{sec:methods-oracle-score}` citation explaining the out-of-band firewall.

2. **Rule 2 — floor-disclosure**. Any `.tex` file containing a claim of the form `"REAL detects X"` where `X` is a canonical rare population name (from catalog v0.3 §2) must either (a) cite `\Cref{sec:minimax-detectability}` and include a `\Cref{sec:discussion-limitations-immunology}` caveat, or (b) be restricted to populations with `recovered_at = flat_whole` in the detectability table.

3. **Rule 3 — evidence-package-gate**. Any `.tex` file claiming a `candidate_novel` population must cite `\Cref{sec:evidence-package}` and provide a companion evidence-card JSON in `workspace/evidence_cards/` that passes the schema.

4. **Rule 4 — hierarchical-transparency**. Any file claiming detection of LAMP3+ mregDC / AS DC / TPEX / KIR+ NK / pDC / FOLR2+ Mφ / TRM-TR / eTreg MUST include the phrase `"hierarchical REAL"` or `"sub-compartment"` in the same paragraph (or in the caption for figures), because these populations are below flat floor per Supp Table S1.

5. **Rule 5 — bib-consistency**. Every `\citep{}` must resolve to an entry in references.bib or in the designated peer bib-additions file (immuno_bib.bib, tumor_bib.bib, etc.).

---

## Version history

- **v0.1 (2026-05-10):** initial cross-reference index covering all immunology LaTeX handoffs as of v0.20 (agent/immunology@07fb8d4).
