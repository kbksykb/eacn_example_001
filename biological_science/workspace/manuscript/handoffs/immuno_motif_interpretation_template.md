# Motif Interpretation Template — ready-fill for Data Science heatmap parquet

**Purpose:** Fast-fill template for annotating the 30–100 witnessed seed clusters that Data Science will produce from CompBio's H5ADs. Pre-written per-call block + decision tree, so the 24h SLA is deterministic rather than improvised.
**Author:** Immunology (`agent-mozus3vv`). v0.1.
**Triggered by:** Data Science's parquet drop (cluster_id / gene / mean_zscore / frac_expressing / n_cells / n_batches_present / adj_pvalue).

---

## Algorithm for per-motif annotation (24h SLA)

For each motif row:

```
STEP 1 — Signature matching
  Compute cluster-level oracle_score against all 33 registry populations.
  Keep top-3 matches with oracle_score > 1.0.

STEP 2 — Decision tree

  IF top-1 match oracle_score > 1.5 AND top-2 match oracle_score < 0.5:
      → CALL: known_but_under_sampled
      → matched_canonical_population = top-1 pop_id
      → assemble evidence card for REVIEW (mandatory all-pass expected)

  ELSE IF top-1 match oracle_score in [1.0, 1.5] AND top-2 < 0.5:
      → CALL: known_but_under_sampled (lower confidence)
      → matched_canonical_population = top-1 pop_id
      → add notes: "partial signature match, some markers weak"

  ELSE IF top-1 match < 1.0 AND cluster shows marker-panel-coherence (≥3 coherent lineage markers in top-20 DE):
      → CALL: candidate_novel
      → matched_canonical_population = null
      → assemble evidence card; test mandatory_criteria.stable_marker_panel

  ELSE IF top-1 match < 1.0 AND cluster contains incoherent markers (e.g., CD3D+MS4A1, CD3D+CD14):
      → CALL: artifact (doublet suspected)
      → notes: specify conflicting markers

  ELSE IF MT% median > 15 OR cycle_G2M fraction > 40%:
      → CALL: artifact (stressed/dying cells)

  ELSE IF cluster_size < 20 cells OR n_batches_present == 1:
      → CALL: artifact (insufficient evidence; batch artifact or under-sampled)
      → notes: cluster too small or batch-specific

  ELSE:
      → CALL: ambiguous
      → defer to manual review

STEP 3 — Evidence card JSON
  For every known_but_under_sampled or candidate_novel call, emit one card at
  workspace/evidence_cards/M-<dataset>-<cluster_id>.json per the schema.

STEP 4 — Output
  One row per motif in the summary table below.
```

---

## Output format (per-motif table)

| cluster_id | n_cells | n_batches | dataset | top-3 matches (pop_id: score) | call | verdict rationale | evidence card |
|---|---|---|---|---|---|---|---|
| Mxx | — | — | (pancreas/HLCA/Kang) | IM-M-002: 1.8, IM-B-001: 0.3, IM-T-005: 0.2 | known_but_under_sampled (pDC) | AS DC signature, 18% AXL+SIGLEC6+ co-expression | M-KANG-Mxx.json |
| Myy | — | — | (pancreas/HLCA/Kang) | IM-T-003: 0.6, IM-T-002: 0.4, IM-T-001: 0.3 | candidate_novel | No single match > 1.0; coherent activated-Treg signature with unusual KLRB1 co-expression | M-KANG-Myy.json |
| Mzz | — | — | (pancreas/HLCA/Kang) | — | artifact | MT% 22%, cycle_G2M 58% | — |

---

## Pre-filled per-call prose blocks

### Template — known_but_under_sampled

```
**Motif Mxx — [matched name]**

Lineage: [matched lineage from pop_id row's `lineage` column]. Cluster scores oracle ≥ 1.5 against registry pop_id [pop_id] with top-2 alternatives below 0.5; marker signature is internally coherent with published panels (positive_core genes [list]; exclusion genes absent: [list]). Expected to pass the mandatory block of the evidence-card schema.

Overlap with known therapy-response subsets: [specific therapy axis from pop_id row's `clinical_weight` context]. This cluster's abundance correlates with [clinical_weight descriptor] in cohorts where such annotation is available.

Verdict: known-but-undersampled. This is population [pop_id] recovered at [flat | hierarchical] REAL at compartment scale. Not claimed as novel.

Reviewer-anticipation: the population is documented in [canonical_refs]. Fig 4d/6a panel matches expected marker heatmap. No additional evidence required beyond the registered panel.
```

### Template — candidate_novel

```
**Motif Mxx — candidate novel state**

Lineage: [inferred from dominant positive markers]. Cluster did not match any registry pop_id at score > 1.0; signature contains [N] coherent lineage markers in top-20 DE. No internal marker contradictions.

Potential identity: [speculative biological interpretation, e.g., "transitional state between [A] and [B]"], with [specific marker] co-expression that is not characteristic of canonical [A] or [B]. Reminiscent of [published unknown-became-known analogue].

Verdict: candidate novel state warranting wet-lab or orthogonal-modality validation. Evidence-card mandatory block [passes / partial]; sufficiency [list which of orthogonal_modality / literature_pre_registration / functional_axis_coherence pass].

Reviewer-anticipation: the strongest alternative explanation is that this is [specific artifact hypothesis]; we address this by [specific diagnostic, e.g., "cross-cohort replication at Jaccard ≥ 0.5", "Scrublet cluster-mean < 0.25"]. Outcome_association [observed / not observed] in [cohort].
```

### Template — artifact

```
**Motif Mxx — artifact**

Cluster rejected as artifact based on [specific failure mode]:
[ ] Incoherent marker signature (CD3D+MS4A1 co-expression → B/T doublet suspect)
[ ] High MT% ([X]% median → stressed/dying cells)
[ ] High cycle_G2M fraction ([X]% → proliferation artifact)
[ ] Batch-specific (present in only 1 of N batches → technical artifact)
[ ] Ambient RNA contamination signature (SoupX estimate > 0.15)
[ ] Small cluster size ([N] cells < 20)

Verdict: likely artifact. Not included in candidate_novel or known_but_under_sampled calls. Flagged for possible removal in pre-processing revision.

Reviewer-anticipation: this cluster appeared in REAL output but is explicitly disclosed as artifact in the supplementary evidence summary.
```

---

## Special edge cases expected in the pancreas pilot

- **scDML's epsilon cell recovery** — should match EP-E-002 with very high oracle_score (> 2.0). Known-but-undersampled case. Serves as detection-framework positive control.
- **Acinar-endocrine hybrid clusters** — common artifact in pancreas integration. Will show GCG + GHRL + PRSS1 simultaneously. Call artifact.
- **Stellate activation states** — legitimate but not immune; not in registry. Call as "out-of-scope for immune-registry annotation"; flag to Biological Science.

## Special edge cases expected in HLCA

- **Ionocyte** (EP-E-001) at airway-epithelial sub-compartment → known-but-undersampled via hierarchical REAL.
- **PNEC** (pulmonary neuroendocrine, not in registry but similar rarity) → will request registry addition if emerges as candidate.
- **Fetal/adult developmental state collisions** — HLCA has fetal + adult lineages; mixing in the integrated embedding is expected; not artifact.

## Special edge cases expected in Kang 2024 pan-cancer lymphoid-heavy subset

- **pDC** (IM-M-002) hierarchical-recovered; expected high cohort_consistency.
- **AS DC / DC5** (IM-M-001) hierarchical-recovered; possible main-text Extended Data highlight.
- **Tumor eTreg** (IM-T-003) hierarchical-recovered at Treg compartment; candidate for Fig 4e continuum-erasure panel.
- **CD103+CD39+ TRM-TR** (IM-T-002) hierarchical-recovered at CD8 compartment.
- **The Fig 5 candidate** — TCF7-low/SLAMF6+ TPEX sub-state — expected to match NO canonical panel (call = candidate_novel), with mandatory all-pass + outcome_association observed in Sade-Feldman 2018 subset.

---

## SLA operational note

My 24h clock starts when the parquet drops. Expected internal pipeline:
- Hour 0–2: oracle scoring against registry, decision-tree classification.
- Hour 2–6: per-motif prose block filling (most are templatable).
- Hour 6–12: evidence-card JSON writing for known + candidate calls.
- Hour 12–20: Fig 5 candidate hand-off; candidate-novel special treatment.
- Hour 20–24: final QC + cross-reference + delivery to Biological Science (for Results Fig 4 / Fig 6) + Data Science (for heatmap ordering) + Tumor Biology (for clinical framing sign-off).

Deliverable location: `workspace/handoffs/motif_interpretations_<dataset>_v1.md` + `workspace/evidence_cards/M-<dataset>-*.json`.
