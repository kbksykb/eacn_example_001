# Experimental Plan v1 — Preserving Previously Unannotated Rare Subpopulations During scRNA-seq Batch Integration

**Author:** Computational Biology Agent (agent-mozur9ik)
**Date:** 2026-05-10
**Status:** draft v1 — integrates inputs from ML (mozurrmd), Data Science (mozutvsy), Biological Science (mozuoaqx), Immunology (mozus3vv), Tumor Biology (mozurh1r), Philosophy (mozur8ub).
**Target venue:** Nature Methods / Nature Biotechnology.

---

## 0. North Star

Produce a **complete CNS-level submission** that satisfies the four published criteria for this problem (Detectability, Protection, Scalability, Computational validation) by the team deadline. Every experiment in this plan exists to fill a specific cell in the figure plan.

## 1. Framework

### 1.1 Hide-and-Measure (Philosophy, mozur8ub)

Evaluation backbone for Detectability and Protection. Protocol:

1. Curated dataset with known rare population(s) (ionocyte, epsilon, MAIT, etc.).
2. Pre-integration: strip rare label; optionally inject synthetic held-out rare clusters of varied geometry.
3. Run method panel; measure detection using label-free / invariant-based signals only.
4. Reintroduce labels ex post to grade the loss rate.

**Constraint (Philosophy):** include at least one dataset whose annotation was NOT produced by a batch-integration pipeline. Baron pancreas (FACS + manual) qualifies; we'll also evaluate against Tabula Sapiens where feasible.

**Negative control (proposed, awaiting Philosophy confirm):** a synthetic "rare mode" that is purely a batch-noise artefact. Framework MUST NOT flag its loss.

### 1.2 Nested 4-Experiment Protocol (Immunology, mozus3vv)

All four experiments are run on each dataset in the panel:

- **E1 Holdout-identity** — real annotated rare subset; hide label, run integration, measure detection.
- **E2 Holdout-cells** — drop X% of rare cells per batch (mimics undersampled rare populations).
- **E3 Dilution series** — inject rare cluster at a sweep of abundance levels (0.01%, 0.05%, 0.1%, 0.5%, 1%).
- **E4 Adversarial mixing** — force batch-specific imbalance of the rare subset (batch A has 5× more than batch B).

Each experiment probes a distinct loss mode.

### 1.3 Detectability Ensemble — "RareScore" (CompBio + ML + TumorBio + DS)

Three complementary channels, combined by late-fusion ensemble per candidate:

1. **Mutual-kNN purity + batch-mixing drift** (CompBio). Pre-integration over-cluster → candidate rare mode C. For each c ∈ C, track: mutual-kNN purity(c), local density(c), per-batch mixing index(c) before and after. Loss if post-integration purity(c) < τ AND NN identity of a majority of c's members is dominated by a different pre-integration cluster.
2. **OT coupling + density-flow matching** (ML, mozurrmd). Conservation of local mass: compute OT plans between pre-integration per-batch manifolds and the post-integration joint embedding; flag cells that lose local mass to abundant neighbours. Permutation null for significance.
3. **Procrustes displacement + bootstrap non-reproducibility** (TumorBio, mozurh1r). Align pre- and post-integration embeddings via Procrustes; measure per-candidate displacement + local neighbourhood-purity increase (i.e., post-integration neighbours pulled from the "wrong" pre-integration cluster); resample batches with replacement — if the candidate is not bootstrap-stable, it is a noise mode, not a rare subpop.

Aggregate metric: **LossRate@k** over top-k rarest candidates, reported with bootstrap 95% CI.

### 1.4 Protection (ML + CompBio)

Two mechanisms, evaluated jointly:

1. **Anchor-preserving integration (CompBio).** Identify rarity-anchors pre-integration (cells with high RareScore candidate membership); during integration, add a contrastive loss keeping each anchor's local neighbourhood stable.
2. **Density-flow regularizer (ML).** Penalize contraction of low-density neighbourhoods; mutual-NN-anchored contrastive alignment so un-anchored rare pockets are not pulled together by batch-alignment.

Drop-in for any integration backbone that exposes an embedding loss (scVI, scDREAMER, scCRAFT).

## 2. Datasets

Staging to `/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/data/`.

### Tier A — Unit tests (small, fast)
- **Pancreas**: Baron, Muraro, Segerstolpe, Wang (~25–40k cells). Rare positive controls: epsilon, activated stellate, Schwann.
  - Baron is manually annotated pre-integration → qualifies for the non-integration-annotated constraint.
- **Lung HLCA core** (Sikkema 2023, ~600k). Rare positive control: pulmonary ionocyte.

### Tier B — Scalability headline
- **Kang 2024 pan-cancer immune atlas** (4.9M cells, 104 publications, 943 patients, 30 cancer types).
  - Pilot: 500k-cell stratified subsample (by cancer type, by batch).
  - Full scale: 4.9M once method panel is stable.

### Tier C — Cancer-specific rare populations (Tumor Biology, mozurh1r)
Evaluate detection on:
- Drug-tolerant persisters, EMT-hybrid, TPEX/TCF1+ CD8, LAMP3+ mregDC, FOLR2+ TAM, CAF subtypes (myCAF / iCAF / apCAF), EC tip cells.

### Tier D — Cross-tissue immunology (Immunology, mozus3vv)
- HCA Immune Cell Atlas + gut atlas — MAIT across tissue batches (Exp-I1).
- Pan-cancer T cell atlas + Kang 2024 — tumor eTreg dilution (Exp-I2).
- COMBAT — pDC severe-COVID prevalence shift (Exp-I3).

### Additional rare anchors (Immunology catalog)
ionocytes, pDC, AS DC/DC5, MAIT, γδ T, TRM CD103+CD39+, eTreg, Tfh/Tfr, ABC, ILCs, TAN, mast, Langerhans, CD34+, epsilon. Oracle marker panels from Immunology, used ex post only.

## 3. Method Panel

Primary: Harmony, Seurat-RPCA, Scanorama, BBKNN, scVI, scANVI, scDML, scDREAMER, scCRAFT. (9 methods.)
Secondary (if bandwidth): CellANOVA, STACAS, RBET.

Each method is evaluated with and without our protection modules plugged in (where applicable).

## 4. Execution Plan

### 4.1 Infrastructure

Shared GPU box `47.103.140.117` (8×A100 80 GB). Layout:

```
/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/
├── data/                                # raw + preprocessed H5AD
├── integrations/<method>/<dataset>.h5ad # X_integrated + obs/var
├── runtime/<method>_<dataset>.json      # wall-clock + peak VRAM + RAM
├── detections/<detector>/<method>_<dataset>.parquet # per-candidate loss calls
└── figures/                             # DS's rendering target
```

Environment: conda-lock pinned (scanpy 1.10, scvi-tools 1.1, harmonypy, seurat-integration-python, scDREAMER, scCRAFT, scDML, Scanorama, BBKNN). Deterministic seeds.

### 4.2 GPU allocation (all 8 A100s in parallel)

- GPU 0–1: pancreas panel (fast).
- GPU 2–3: HLCA ionocyte panel.
- GPU 4–5: pan-cancer 500k pilot (scVI / scANVI / scDREAMER training).
- GPU 6: Tumor-Biology cancer-rare evaluations.
- GPU 7: ML's density-flow detector training + RareScore scoring.

Policy (from repo CLAUDE.md): if a GPU has free VRAM ≥ a queued job's working set, launch now. Don't reserve headroom.

### 4.3 Phased schedule

- **Phase 1 (this cycle, pancreas pilot):** download pancreas datasets, run 9-method panel on Baron+Muraro+Segerstolpe+Wang, apply RareScore ensemble on epsilon holdout, confirm detection (Harmony/Seurat should flag, scVI/scDREAMER should not).
- **Phase 2:** HLCA ionocyte, same protocol.
- **Phase 3:** synthetic sweep (Experiments E2/E3/E4 on pancreas), factorial knobs (abundance, mixing, geometry, UMI depth).
- **Phase 4:** pan-cancer 500k pilot; cancer-specific rare evaluations with Tumor Biology.
- **Phase 5:** protection mechanism — add anchor-preserving + density-flow regularizer to scVI/scDREAMER/scCRAFT; re-run all tiers, show LossRate@k drop without scIB-E major-type degradation.
- **Phase 6:** full 4.9M Kang 2024 — scalability demonstration.

## 5. Deliverables

- `workspace/code/rarescore/` — detector library (CompBio + ML + TumorBio channels, ensemble).
- `workspace/code/harness/` — unified integration runner (DS-compatible LRS output).
- `workspace/code/protection/` — anchor-preserving + density-flow loss modules.
- `workspace/results/` — per-phase parquet + YAML configs.
- `workspace/figures/` — paper figures (DS owns rendering).
- `workspace/manuscript/` — LaTeX (BioScience owns final integration).

## 6. Open Questions (awaiting peer input)

- Philosophy (mozur8ub): confirm the batch-artefact negative control.
- ML (mozurrmd): concrete reference papers for the density-flow regularizer.
- Biological Science (mozuoaqx): rigor/manuscript lens on the plan — any reshape before GPU burn.
- Mathematics: no agent yet registered. Will need for LossRate@k theoretical bounds, exchangeability assumptions, and any concentration results for Procrustes-bootstrap stability. CompBio will seed a create_task once Mathematics agent comes online.

## 7. Risk register

| Risk | Mitigation |
|---|---|
| Kang 2024 download / licensing delay | Start with CELLxGENE mirror; parallel subsampling from HLCA if unavailable. |
| RareScore channels disagree | Pre-register combiner (Fisher / Stouffer) on pancreas before scaling. |
| Protection hurts major-cell-type batch correction | scIB-E suite is run on every protected run; treat as a hard guard. |
| GPU contention with DS/BioScience | All three are GPU-priority; use `CUDA_VISIBLE_DEVICES` per CLAUDE.md and co-locate small jobs. |
