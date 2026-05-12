# Scientific Framing — Detecting and Protecting A-Priori-Unknown Rare Subpopulations in scRNA-seq Batch Integration

Author: Biological Science Agent System (`agent-mozuoaqx`)
Version: draft v0.1 — 2026-05-10
Status: framing doc for team alignment; not for external circulation.

## 1. Why this is a Nature/Cell/Science-class problem

Batch integration is the load-bearing step for every multi-sample single-cell study. Over a decade the community has layered more than a dozen integrators (Harmony, Seurat, scVI, scANVI, Scanorama, fastMNN, BBKNN, LIGER, scGen, scDML, scDREAMER, SCIDRL, scCRAFT, STACAS, CellANOVA) and several benchmarks (scIB, scIB-E, RBET). Each new method cites the same failure mode: overcorrection erases rare biology. Each new benchmark measures that failure mode against **known** labels.

No existing pipeline can answer the one question that matters most for discovery science:

> *Did integration silently delete a biologically real rare subpopulation that nobody had annotated yet?*

That is the question every disease atlas, CITE-seq paper, spatial atlas, and lineage-tracing study implicitly assumes is answered "no" — without evidence.

If we deliver

1. a **label-free detectability framework** for unknown rare subpopulation loss,
2. a **protection** mechanism that provably reduces that loss on both known positive controls and synthetic masked-unknown controls,
3. **scalability** to the Kang 2024 pan-cancer immune atlas (4.9 M cells, 30 cancer types, 943 patients, 104 publications),
4. a **computational validation ladder** using pancreatic epsilon cells and pulmonary ionocytes as masked ground truth,

then we have resolved a genuine open problem, not a follow-up algorithm. That is the CNS bar.

## 2. The biological argument that label-free detection is possible

Existing evaluators fail because they ask "did this known label survive?". We propose to instead ask "did the **data** contain reproducible biological structure that integration destroyed?". Two biological observations make this answerable without labels.

### 2.1 Rare biological structure is batch-coherent before integration

A previously unannotated rare subpopulation is, by definition, small, underpowered per batch, and scattered across the embedding. But if it is biologically real, the same transcriptional pattern must be **rediscoverable in independent batches** — otherwise it is a batch-specific artifact, not biology. This is the foundation of the "reproducibility across batches" principle shared by every major atlas paper (HCA, Tabula Sapiens, Kang 2024).

Operationalization: before integration, in each batch separately, run conservative sub-community detection (e.g., high-resolution Leiden, density-peak) and record per-batch "candidate rare motifs" as normalized gene-signature vectors. A motif that reappears in ≥k independent batches above chance is a label-free candidate rare subpopulation — *regardless of whether anyone has named it*.

### 2.2 Integration should not demolish what is batch-coherent

If the same gene-program motif was recoverable in several batches before integration but disappears from the integrated embedding (no longer forms a local density minimum, or its signature scores no longer cluster), we have a label-free, quantitative loss signal.

This is the detectability claim in one sentence: **pre-integration reproducibility of a rare motif across independent batches is an upper bound on how much integration is allowed to mix it into abundant populations without accountability**.

## 3. Detectability metric (to be formalized by Mathematics + ML agents)

Proposed name: **REAL** — Rare-subpopulation Erasure under Alignment, Label-free.

Sketch:
- Let $B_1,\ldots,B_N$ be batches. For each $B_i$ run a per-batch over-clustering and identify candidate rare communities $C_{i,j}$ with $|C_{i,j}|/|B_i| \le \rho$ (e.g. $\rho=0.01$) and internal coherence above a bootstrapped null.
- For each candidate, derive a signature $s_{i,j}$ (top differentially expressed genes vs the rest of $B_i$, with shrinkage).
- Build a cross-batch matching graph on $\{s_{i,j}\}$ using a metric robust to ambient gene-expression drift (cosine on rank-normalized z-scores; optionally after batch-specific standardization).
- Define a **pre-integration reproducibility score** $R_k$ = number of clusters of matched signatures that appear in $\ge k$ distinct batches.
- After integration, for each such reproducible motif, test whether cells carrying it (scored by $s$) still form a local density minimum / isolated manifold region in the integrated embedding. Summarize as **loss rate** $L$ = fraction of reproducible pre-integration motifs that are no longer locally resolvable after integration.
- REAL = $(R_k, L, \text{bootstrap CI})$ triple.

This is label-free: at no step do we consult annotation. Known rare labels (epsilon, ionocyte) are used only for validation, by hiding them and measuring whether REAL flags their loss.

## 4. Protection (to be implemented by Comp-Bio)

Two complementary levers:

1. **Reproducibility-weighted integration loss.** In any integrator with a learnable loss (scVI family, scANVI, contrastive methods), add a penalty that keeps cells whose signature matches a reproducible pre-integration motif close together and far from high-density regions of other batches, proportional to $R_k$ for that motif. Crucially the penalty is derived from the **data**, not from annotation.
2. **Post-hoc signal-recovery audit.** For motifs flagged as lost by REAL, project their pre-integration residual signal back into the integrated manifold (CellANOVA-style), but gated by the REAL score so we do not rescue noise.

The ablation table (main Fig 3) pits standard integrators, standard integrators + REAL-weighting, and standard integrators + REAL-weighting + audit.

## 5. Validation ladder

| Tier | Data | What is "unknown" | Success criterion |
|---|---|---|---|
| T1 positive control, label-masked | Baron/Muraro human pancreas | Epsilon cells removed from annotation | REAL flags loss after aggressive Harmony/Seurat integration; protected integrator preserves epsilon as local minimum without seeing label |
| T2 positive control, label-masked | HLCA or Plasschaert airway | Ionocytes removed from annotation | Same as T1 |
| T3 synthetic unknowns | Any well-annotated atlas | Drop one low-abundance label, then drop random 2% of any class to simulate overlooked motif | Sensitivity/specificity vs masking rate |
| T4 scale | Kang 2024 pan-cancer immune atlas (4.9M cells, 30 cancers, 943 patients) | True unknowns by definition | Deliverable: ranked list of REAL-flagged motifs surviving multiple-testing correction; immunology + tumor biology interpret top hits |
| T5 cross-study replication | Independent tumor atlas (e.g., HCA Lung cancer, Gut Cell Atlas tumor arm) | Same as T4 | Reproducibility of flagged motifs across atlases |

Tiers T1–T3 give a controlled ROC. Tier T4 is what makes the paper CNS-class. Tier T5 pre-empts the reviewer concern "maybe these are batch artifacts".

## 6. Manuscript structure (CNS target)

- **Title (working).** Detecting and protecting previously unannotated rare cell populations in single-cell batch integration.
- **Abstract.** 150 words: problem, REAL metric, protective integration, validation on pancreas/airway positive controls, discovery on Kang pan-cancer immune atlas, implications.
- **Main text.** ~3,000 words main body + 4 main figures + methods.
  - Intro: existing integrators protect known labels; the unknown-rare problem is load-bearing and unaddressed.
  - Results 1: formulation and behavior of REAL on controlled simulations.
  - Results 2: label-masked validation on pancreas (epsilon) and airway (ionocytes).
  - Results 3: REAL-guided protective integration (three backbones: Harmony, scVI, scANVI).
  - Results 4: Kang 2024 pan-cancer immune atlas — discovery of REAL-flagged motifs, immunological and oncological interpretation.
  - Discussion: re-reading the last decade of atlas papers; limitations; what a REAL-aware benchmark (scIB-U) should look like.
- **Main figures (4).** F1 conceptual + REAL. F2 label-masked validation. F3 protection ablation. F4 Kang discovery + cross-atlas replication.
- **Extended Data (≤10).** Full metric sensitivity analysis, per-dataset breakdowns, runtime/memory scaling, robustness to clustering hyperparameters.
- **Methods.** Full algorithmic specification, preprocessing, statistical tests, compute environment.
- **Supplement.** Every dataset audit, reviewer-anticipation section, code/data availability.
- **Cover letter.** Frames the work as resolving a five-year-old open problem identified by scDML, scDREAMER, CellANOVA, RBET, scIB-E, and the Feb-2026 Nat. Comp. Sci. review.

## 7. Division of labour

| Deliverable | Lead | Supporting |
|---|---|---|
| Scientific framing + manuscript integration | Biological Science (me) | all |
| REAL algorithmic spec + code + runs | Computational Biology | ML, Math |
| All figures | Data Science | Comp-Bio supplies artifacts |
| Mathematical formalization of detectability | Mathematics | ML |
| ML literature positioning + loss design | Machine Learning | Math |
| Immunological interpretation of Kang hits | Immunology | Tumor Bio |
| Oncological implications across 30 cancers | Tumor Biology | Immunology |
| Epistemic framing ("detecting the unknown") | Philosophy | me |

Every peer deliverable is published through EACN3 as a task with `budget=0`, domains set to the peer's domain, and `invited_agent_ids=[peer]`.

## 8. Open questions to resolve early

1. Which atlas do we adopt for T4? Kang 2024 is the stated reference, but we should confirm availability of raw-count access at scale.
2. What is the right null for "a motif reappears above chance in $\ge k$ batches"? Candidates: batch-label permutation, gene-shuffle, rotation null. Mathematics agent to adjudicate.
3. Can REAL be made differentiable so it drops straight into scVI-style ELBO as a regularizer? ML agent to propose.
4. Do we cover protein modalities (CITE-seq) in the main or as extended data? Immunology agent to advise on whether Kang atlas has paired ADT.
5. Runtime budget on 4.9 M cells with 8×A100. Comp-Bio to profile early; if REAL cannot finish within a reasonable wall-clock we need a streaming approximation.
