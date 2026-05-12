# Computational Biology Agent — State (2026-05-12)

## Agent Identity
- agent_id: agent-mozur9ik
- Network: http://127.0.0.1:8888 (connected, claimed)
- Loop: CronJob 2acf174a, every 20 min

## Criteria Status
- **Criterion 1 (Detectability)**: DECISIVELY MET
- **Criterion 2 (Protection)**: DECISIVELY MET
- **Criterion 3 (Scalability)**: FULLY MET (this session)
- **Criterion 4 (Computational validation)**: MET

## Kang 2024 Fig 4 — COMPLETE (agent/computational_biology@e46afa3)

### 5-Compartment Harmony Detection Summary
| Motif ID | Compartment | Label | n_cells | pi | Δ | tau | p | sig |
|---|---|---|---|---|---|---|---|---|
| M-KANG-MYL-001 | myeloid | LAMP3+ mregDC | 225 | 0.010 | 3.337σ | 93.82 | 0.005 | ✓ |
| M-KANG-TNK-001 | T/NK | TPEX (canonical, Miller 2019) | 314 | 0.020 | 1.903σ | 84.61 | 0.005 | ✓ |
| M-KANG-B-001 | B | plasmablast | 187 | 0.020 | — | 64.32 | 0.005 | ✓ |
| M-KANG-MES-001 | mesenchymal | myCAF | 447 | 0.020 | — | 63.16 | 0.005 | ✓ |
| M-KANG-EPI-001 | epithelial | ionocyte | 291 | 0.010 | — | 34.86 | 1.000 | ✗ |

4/5 detected. Ionocyte negative = valid specificity control.

**IM.4 compliance note (Philosophy flag)**: M-KANG-TNK-001 is **hierarchical-recovered** detection within the T/NK sub-compartment — NOT above the flat whole-atlas floor. Any manuscript claim of TPEX detection must name "hierarchical REAL / T/NK sub-compartment" context. Flat whole-atlas REAL does not attain this detection.

**Biological identity clarity (Philosophy flag)**: Kang T/NK markers TCF7/TOX/PDCD1/CXCR5 = **canonical TPEX (Miller 2019)**. This is DISTINCT from the Zheng MM CXCL13+ terminal Tex result (prior session). §05/§07 manuscript text must make this distinction explicit.

### Multi-Method on Kang Myeloid (method-agnostic detection)
- Harmony: p=0.005, τ=93.82
- scVI:    p=0.005, τ=124.51
- Scanorama: incomplete (server restart)

### Cancer-Type Specific (myeloid, single CancerAbbr)
- PAAD (n=1,076): p=0.005, τ=138.07 (strongest tau)
- CRC (n=1,383): p=0.015, τ=94.60
- HNSC (n=1,883): p=0.149, τ=74.22 (batch-limited)

### Cross-Atlas Replication (Salcher 2022 NSCLC)
- Full atlas (66,571 cells, 36 batches): p=0.97, τ=4.57 — NOT detected
- Myeloid-only (13,315 cells, 36 batches): p=0.005, τ=9.16 — DETECTED
- **Amplification 5.0× validates hierarchical REAL protocol**

### Three-Atlas LAMP3+ mregDC Confirmation
1. Kang 2024 pan-cancer myeloid: p=0.005
2. Cheng 2021 PAAD myeloid: p=0.01
3. Salcher 2022 NSCLC myeloid: p=0.005

### Key Theoretical Validation (for Math Theorem S1-hierarchical)
- Myl: Δ=3.337σ, π=0.01, n·π²·Δ²=25.08 (cluster-regime, hierarchical-recovered at myeloid sub-compartment)
- TNK: Δ=1.903σ, π=0.02, n·π²·Δ²=22.75 (state-regime, hierarchical-recovered at T/NK sub-compartment)
- Both above the **compartment-level** floor; neither is above the flat whole-atlas floor
- Detection requires hierarchical REAL (compartment-stratified) — flat whole-atlas REAL would not attain these
- Kang amplification 4.4× (myl), 6.3× (tnk)
- Salcher 5.0× — cleanest demonstration (same dataset: flat fails p=0.97, compartment-stratified fires p=0.005)

### τ-Scaling (workspace/results/kang2024_fig4/tau_scaling.json)
- B∈{100,300,500,750,1043}: all hit permutation floor (p=0.005)
- τ varies 88-110 — robustness across scales confirmed
- For √B scaling validation: use Panel A synthetic sweep (π≈0.005, Δ≈0.8σ)

### D-Column Per-Cancer-Type (myl LAMP3+)
SSCC 9.05% (+0.081), HNSC 5.84% (+0.048), CRC 2.89% (+0.019) — enriched
PAAD 0.28% (-0.007), LC 0.38% (-0.006) — depleted

### Melanoma Negative Results (Trap S0.8 — envelope framing per Philosophy flag)
- mel_GSE200218 cutaneous (28k cells, 5 batches): TPEX p=0.149, LAMP3 p=0.388
- Merged 4-cohort melanoma (144k cells, ≤20 batches): TPEX p=1.0, LAMP3 p=0.46

**Correct framing for manuscript** (per Trap S0.8): these are valid negative results meaning **"no motif above detection threshold at this scale"**, NOT "these populations are absent." At n=28k / 5 batches the Theorem 1 minimax envelope (n·π²·Δ² ≥ floor) may not be met — the detection floor rises with fewer batches, and the signal is below envelope. The same populations are detected at pan-cancer scale (Kang T/NK TPEX p=0.005 at 953 batches; Kang myeloid LAMP3+ p=0.005 at 1,043 batches).

## Git Clone
- ~/eacn_cb_clone (use this, NOT ~/eacn_example_001)
- .gitignore: /biological_science/, /data/, /results/, *.npy, *.h5ad

## Pending (Next Cycle)
1. Scanorama on Kang myeloid (conservative, single-job, low-load)
2. 15-study atlas merge (conservative, low-parallelism)
3. Baron/Muraro/Segerstolpe/Wang pancreas (§08 un-caveat — blocked on data access)
4. Caushi 2021 NSCLC (separate GEO download)
5. Panel A tau-scaling row selection (coordinate with DataSci)

## CPU/Memory Policy (per user instruction)
- Max 1-2 CPU-heavy jobs at a time
- OT permutations on subsampled n<20k
- Skip multi-batch Scanorama on >500 batches
- Monitor `uptime` before launching; wait if load > 8


5. Caushi 2021 NSCLC download (GSE173351)

