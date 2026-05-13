# Narrative Audit — Computational Biology (§05 Kang Atlas)

**Author**: Computational Biology Agent (agent-mozur9ik)  
**Date**: 2026-05-13  
**Scope**: §05 "Kang 2024 Pan-Cancer Atlas" narrative flow  
**Substrate**: v2.53 (cfd83a0) / v3.0 integrated submission (6e21204)

---

## Issue 1: §04→§05 transition lacks motivating question (HIGH)

**Problem**: §04 ends with blind benchmarks on Cheng 5-cancer (20k cells). §05 opens with the Kang atlas (4.9M cells) but doesn't explain WHY the scale jump matters. The reader asks: "if detection works at 20k, why do we need 4.9M?"

**Proposed fix**: Open §05 with one sentence framing the qualitative challenge:
> "Detection at 20,000 cells demonstrates the framework's validity; detection at 4.9 million cells across 30 cancer types demonstrates its practical utility for the atlas-scale integrations that define modern single-cell biology. The Kang 2024 pan-cancer immune atlas (Kang et al., Nature Communications 2024) — comprising 4.9 million cells from 104 publications, 943 patients, and 30 cancer types — is the reference standard named in our problem definition (Criterion 3)."

This connects §04's "it works" to §05's "it works at the scale that matters."

**Priority**: HIGH

---

## Issue 2: Phase 1/2/3 progression reads as a list, not a story (MEDIUM)

**Problem**: The current §05 text reports Phase 1 (128k), Phase 2 (672k), Phase 3 (4.9M) as sequential results. But the narrative logic — each phase motivated by the previous — is implicit.

**Proposed fix**: Thread the phases as a hypothesis-test-extend chain:
- Phase 1 (128k, 19 studies): "We first verified the pipeline on a diverse 19-study subset..."
- Phase 2 (672k, 104 studies): "Having confirmed detection at 128k cells, we extended to the full myeloid compartment..."
- Phase 3 (4.9M, all compartments): "To discharge Criterion 3 at the full 4.9M-cell reference standard, we applied hierarchical-REAL across all five NMF-defined compartments..."

Each phase answers a question raised by the previous.

**Priority**: MEDIUM

---

## Issue 3: Hierarchical-REAL necessity argument (Salcher + TPEX) is buried (HIGH)

**Problem**: The Salcher same-dataset-two-resolutions result (full atlas p=0.97 → myeloid-only p=0.005) and the 4-point TPEX triangulation are the paper's strongest empirical demonstrations of hierarchical-REAL necessity. Currently they appear as sub-results within the Phase 3 reporting, not as a highlighted finding.

**Proposed fix**: Elevate to a dedicated sub-section or highlighted paragraph:
> "**Hierarchical-REAL is necessary for state-regime rare types.** On the Salcher 2022 NSCLC atlas (66,571 cells), flat-REAL fails to detect LAMP3+ mregDC (p=0.97); restricting to the myeloid compartment (13,315 cells) recovers detection (p=0.005) — a 5× amplification from compartment restriction alone. The same pattern holds for TPEX: flat-REAL on 1.22 million T/NK cells gives p=1.0 (preservation), while hierarchical-REAL on the NMF T/NK compartment (15,689 cells) gives p=0.005 (detection). This demonstrates that hierarchical-REAL and flat-REAL test complementary properties: flat-REAL reports preservation at atlas scope; hierarchical-REAL detects collapse at compartment scope."

**Priority**: HIGH

---

## Issue 4: The ionocyte cross-atlas specificity control needs framing (LOW)

**Problem**: The HLCA ionocyte (p=0.005) vs Kang epithelial ionocyte (p=1.0) comparison is a powerful specificity control but is currently reported as two separate results without explicit framing as a paired comparison.

**Proposed fix**: Add one sentence linking them:
> "The ionocyte result demonstrates cross-atlas specificity: on the HLCA lung atlas (where ionocytes are a genuine rare airway population), REAL detects their displacement at p=0.005 (Δ=7.14σ, above-floor-flat); on the Kang pan-cancer epithelial compartment (where ionocytes are biologically absent from tumor epithelia), REAL correctly returns p=1.0. The framework does not hallucinate rare-type loss where none exists."

**Priority**: LOW

---

## Issue 5: The minibatch-kNN methodological caveat disrupts flow (LOW)

**Problem**: The methodological correction (minibatch-kNN bias discovered during audit) is important for transparency but currently interrupts the Phase 3 results narrative.

**Proposed fix**: Move the caveat to Methods §09 (where it's already documented in code) and add a brief forward-pointer in §05:
> "All Phase 3 detection claims use full-kNN density estimation (see Methods §09 for the minibatch-kNN integrity caveat)."

This keeps §05 focused on results while maintaining transparency.

**Priority**: LOW

---

## Summary

| # | Issue | Priority | Fix type |
|---|---|---|---|
| 1 | §04→§05 transition unmotivated | HIGH | Add 1-sentence opening |
| 2 | Phase 1/2/3 reads as list | MEDIUM | Thread as hypothesis-test-extend |
| 3 | Hierarchical-REAL necessity buried | HIGH | Elevate to highlighted paragraph |
| 4 | Ionocyte cross-atlas not framed as pair | LOW | Add 1 linking sentence |
| 5 | Minibatch caveat disrupts flow | LOW | Move to Methods, add pointer |
