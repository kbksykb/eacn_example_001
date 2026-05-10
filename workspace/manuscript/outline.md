# Manuscript Outline — Detecting and Protecting Previously-Unannotated Rare Subpopulations in Single-Cell Batch Integration

**Working title:** "Label-free detection and preservation of unannotated rare subpopulations in single-cell batch integration."
**Target venue:** Nature Methods (primary) / Nature Biotechnology (secondary).
**Draft holder:** Biological Science agent (agent-mozuoaqx) integrates; CompBio (this agent) supplies methods + results.

---

## Title (working options)

1. "Label-free detection and protection of previously-unannotated rare subpopulations in single-cell batch integration."
2. "When integration hides what you never knew you had: a detectability framework for rare cell-type loss."

## Abstract (150–180 words, placeholder)

Batch integration is universal in multi-sample single-cell analysis and known to overcorrect, erasing rare subpopulations. Every published integration benchmark evaluates preservation of *annotated* rare types; no framework detects loss of rare subpopulations that are unannotated a priori. We present **RareScore**, a label-free detectability framework combining OT mass-conservation, mutual-kNN purity, and Procrustes-aligned displacement signatures stabilized by bootstrap, together with a drop-in **protection module** (anchor-preserving contrastive loss + density-flow regularizer) that attaches to modern integration backbones. On pancreas (epsilon) and lung (ionocyte) hold-out benchmarks, RareScore recovers the known-vulnerable integrator behaviours with AUC > 0.9 and without any external label. Protection reduces rare-subpopulation loss rate by [TBD]% on the 4.9M-cell Kang 2024 pan-cancer immune atlas while preserving scIB-E major-type correction within 5% of unprotected baselines. The framework reveals [N] previously-unreported rare TME subpopulations lost by standard pan-cancer integrations.

## Introduction (Figure 1)

1. **Batch integration is foundational** — 10+ years of methodology, scIB establishing benchmarks, every major atlas uses it.
2. **Overcorrection erases rare subpopulations** — cite scDML, scDREAMER, SCIDRL, scCRAFT, STACAS, CellANOVA, RBET, scIB-E, Feb 2026 Nat Comp Sci review.
3. **But every existing evaluation needs labels.** scDML validates epsilon because epsilon is already annotated. No metric exists for loss of rare types that were never annotated. (Use Philosophy's "latently discoverable" operational framing here.)
4. **Consequence:** all downstream discoveries from integrated atlases rest on an untested assumption.
5. **Our contribution:** (a) label-free detectability framework RareScore; (b) drop-in protection module; (c) demonstration on pancreas, lung, and a 4.9M-cell pan-cancer atlas.

**Figure 1.** Schematic: annotated vs unannotated rare subpopulations; how existing metrics become blind; our proposed detection + protection loop.

## Results

### Result 1 — RareScore detects rare-subpopulation loss without labels (Figure 2)

Pancreas epsilon and lung ionocyte hold-out. 9 integrators. RareScore AUC for calling known-vulnerable vs known-preserving methods. Ablation across the three channels (OT, mutual-kNN, Procrustes).

**Figure 2.**
- (a) Protocol schematic (hide-and-measure).
- (b) Per-method LossRate@k on pancreas (9 bars, epsilon as ground truth).
- (c) Per-method LossRate@k on lung (ionocyte).
- (d) ROC of each RareScore channel for known-loss calls.
- (e) Ensemble vs single-channel improvement.

### Result 2 — Synthetic factorial sweep characterizes the framework's operating envelope (Figure 3)

Knobs: abundance (0.01% → 5%), batch-mixing ratio, geometric proximity to the nearest abundant population, per-cell UMI depth. Heatmap of detection power vs each knob.

**Figure 3.**
- (a) Abundance-proximity detection heatmap.
- (b) Depth-mixing detection heatmap.
- (c) Adversarial-mixing experiment (E4): detection under batch-specific imbalance.
- (d) Negative control: purely batch-artefactual mode is correctly not flagged.

### Result 3 — Protection module reduces loss rate without hurting major-type correction (Figure 4)

Apply anchor-preserving + density-flow to scVI / scDREAMER / scCRAFT. Evaluate with RareScore and full scIB-E.

**Figure 4.**
- (a) Paired before/after LossRate@k per method.
- (b) scIB-E major-type batch correction (should be unchanged or slightly improved).
- (c) UMAP insets showing preserved rare clusters.

### Result 4 — 4.9M-cell pan-cancer scalability (Figure 5)

Kang 2024. Runtime, peak VRAM, accuracy vs 500k pilot. Detected vs preserved rare populations in TME.

**Figure 5.**
- (a) Scaling curve (cells vs wallclock), multi-GPU.
- (b) Peak-VRAM vs dataset size per method.
- (c) UMAP of the 4.9M atlas with RareScore-flagged rare pockets highlighted.
- (d) Table: detected rare TME populations before vs after protection.

### Result 5 — Biological discovery (Figure 6)

Cancer-specific rare populations detected under protection that were lost without it. Marker enrichment, cross-tissue consistency, clinical correlate (if available). Tumor Biology + Immunology drive the interpretation.

**Figure 6.**
- (a) Marker expression of newly-preserved rare TME subpopulations.
- (b) Cross-cancer prevalence.
- (c) Survival / response correlate for at least one newly-preserved subset.

## Discussion

- Operational redefinition of "unannotated rare subpopulation" as "latently discoverable in pre-integration signal" — Philosophy's framing.
- Comparison to scIB-E and RBET: RareScore detects a signal orthogonal to both.
- Limitations: framework requires a usable pre-integration embedding; datasets with extreme technical dominance (e.g., scATAC-only) require adaptation.
- Implications for every published atlas that used a non-protected integrator: how many unannotated rare populations were silently lost?

## Methods

- Datasets (preprocessing, batch definitions, QC).
- Integration method configurations (exact versions, seeds).
- RareScore formal definitions (OT, mutual-kNN, Procrustes, bootstrap).
- Protection module derivation (loss gradients, scaling).
- scIB-E evaluation setup.
- Compute (8×A100 box, 80 GB per card).
- Reproducibility (conda-lock, seeds, YAML).

## Supplement

- Per-method per-dataset detailed metrics.
- Factorial sweep full grid.
- Synthetic injection protocol.
- Oracle marker-panel scores vs RareScore calls (validation only, never fed into integration).
- Runtime + VRAM tables for every run.
- Ablations on each RareScore channel and on the protection loss components.

## Cover letter (BioScience drafts)

Key claim: **first framework that solves the criteria in the problem statement — label-free detection + protection + scale to 4.9M.**

## Author list

8 discipline agents + human PI. Equal-contribution for top 3 GPU-primary agents (CompBio, DS, BioScience); corresponding via BioScience (manuscript integrator).
