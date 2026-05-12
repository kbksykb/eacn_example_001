# Computational Biology Agent — State (2026-05-12)

## Agent Identity
- agent_id: agent-mozur9ik
- Network: http://127.0.0.1:8888 (connected, claimed)
- Loop: CronJob 2acf174a, every 20 min

## Prior Work Summary (v3.0-rc, commit 471a2f0)

### Criteria Status
- **Criterion 1 (Detectability)**: DECISIVELY MET
  - Synth factorial, retina abundance series, Cheng PAAD, Cheng5 pan-cancer
  - Panel A v2: 18 operating points, Theorem 1 validated
- **Criterion 2 (Protection)**: DECISIVELY MET
  - RareShield A/B on synth, retina, Cheng PAAD, Cheng5
  - Ablation grid: full +0.361 ΔAUPRC, −L_mass 0.000, F=identity −0.079
- **Criterion 3 (Scalability)**: PARTIALLY MET → IN PROGRESS (Kang 2024)
  - Cheng5 20k cells done; Kang 2024 4.9M now unblocked
- **Criterion 4 (Computational validation)**: MET
  - Pancreas epsilon (synth proxy), retina BC subtypes, LAMP3+ mregDC

## Current Session Work (2026-05-12)

### Kang 2024 Fig 4 Pipeline
- Data extracted: `/mnt/kang_2024_pan_cancer_atlas/extracted/`
  - NMF h5ads decompressed: b, myl, tnk, mesenchymal (epi missing from tar)
  - Atlas dataset: ~33+ per-study h5ads (extraction ongoing)
- Pipeline script: `workspace/code/pilots/kang2024_fig4_pipeline.py`
- Launch script: `workspace/code/pilots/launch_kang2024_fig4.sh`
- Results dir: `workspace/results/kang2024_fig4/`

### Compartment Run Status
| Compartment | n_obs | n_batches | t_total | Status |
|-------------|------:|----------:|--------:|--------|
| b (B cells) | 9,334 | 662 | 22.4s | ✓ DONE |
| tnk (T/NK) | 15,689 | 953 | 47.9s | ✓ DONE |
| mesenchymal | 22,310 | ? | 76.8s | ✓ DONE |
| myl (myeloid) | 22,472 | ? | IN PROGRESS | OT detection running |
| epi (epithelial) | ? | ? | NOT STARTED | .xz not in tar yet |

### Key Finding So Far
- Harmony integration scales cleanly: 9k→22k cells, 662→953 batches, <80s per compartment
- Batch key: `Patient_Organ_Tissue` (consistent across all compartments)
- LAMP3+ detection running on myl compartment (GPU 1, 6.3GB VRAM)

## Remaining Blockers
1. Myl OT detection completing (in progress)
2. Epi compartment: need to check if epi_NMF.h5ad.xz is in the full tar
3. Full atlas h5ads: need to assess total cell count across all per-study files
4. Fig 4 summary report: compile scalability table once all compartments done

## Next Steps
1. Wait for myl to complete → collect detection result
2. Check atlas_dataset extraction for total cell count
3. Write Fig 4 scalability summary (workspace/results/kang2024_fig4_scalability.md)
4. Notify team via EACN3 with results
5. Commit to agent/computational_biology branch
