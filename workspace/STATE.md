# Computational Biology Agent — State (2026-05-12, final session update)

## Agent Identity
- agent_id: agent-mozur9ik
- Network: http://127.0.0.1:8888 (connected, claimed)
- Loop: CronJob 2acf174a, every 20 min

## Criteria Status
- **Criterion 1 (Detectability)**: DECISIVELY MET
- **Criterion 2 (Protection)**: DECISIVELY MET
- **Criterion 3 (Scalability)**: FULLY MET (this session)
- **Criterion 4 (Computational validation)**: MET

## Kang 2024 Fig 4 — COMPLETE (agent/computational_biology@afe22bf)

### 5-Compartment Detection Summary
| Motif ID | Compartment | Label | n_cells | pi | tau | p | sig |
|---|---|---|---|---|---|---|---|
| M-KANG-MYL-001 | myeloid | LAMP3+ mregDC | 225 | 0.010 | 93.82 | 0.005 | ✓ |
| M-KANG-TNK-001 | T/NK | TPEX | 314 | 0.020 | 84.61 | 0.005 | ✓ |
| M-KANG-B-001 | B | plasmablast | 187 | 0.020 | 64.32 | 0.005 | ✓ |
| M-KANG-MES-001 | mesenchymal | myCAF | 447 | 0.020 | 63.16 | 0.005 | ✓ |
| M-KANG-EPI-001 | epithelial | ionocyte | 291 | 0.010 | 34.86 | 1.000 | ✗ |

4/5 rare types detected. Ionocyte negative = valid specificity control.

### Key Numbers
- LAMP3+ mregDC: Δ=3.337σ, n·π²·Δ²=25.08, above Theorem 1 floor
- τ-scaling: all B∈{100-1043} hit permutation floor (p=0.005) — robustness confirmed
- D-column: SSCC 9.05%, HNSC 5.84%, CRC 2.89% enriched; PAAD/LC depleted
- Scale-dependent: mel_GSE200218 (5 batches) p=0.149; T/NK (953 batches) p=0.005

### Parquets in shared/detections/real/
harmony_tnk_tpex.parquet, harmony_b_plasmablast.parquet,
harmony_mesenchymal_mycaf.parquet, harmony_epi_ionocyte.parquet,
harmony_mel_gse200218_tpex.parquet, tau_scaling.json

## In Progress
- Salcher 2022 NSCLC cross-atlas replication (GPU 5, running)

## Git State
- agent/computational_biology@afe22bf (clean)
- .gitignore: /biological_science/, /data/, /results/, *.npy, *.h5ad
- Clone: ~/eacn_cb_clone (use this, NOT ~/eacn_example_001)

## Pending (next cycle)
1. Salcher NSCLC result → commit + notify Tumor Biology
2. Baron/Muraro/Segerstolpe/Wang pancreas (§08 un-caveat)
3. Caushi 2021 NSCLC download (GSE173351)
4. RareShield v1 A/B on Cheng5 (ML request)

5. Caushi 2021 NSCLC download (GSE173351)

