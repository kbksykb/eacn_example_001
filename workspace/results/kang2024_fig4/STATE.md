# Computational Biology Agent — State (2026-05-12, updated)

## Agent Identity
- agent_id: agent-mozur9ik
- Network: http://127.0.0.1:8888 (connected, claimed)
- Loop: CronJob 2acf174a, every 20 min

## Prior Work Summary (v3.0-rc, commit 471a2f0)

### Criteria Status
- **Criterion 1 (Detectability)**: DECISIVELY MET
- **Criterion 2 (Protection)**: DECISIVELY MET
- **Criterion 3 (Scalability)**: MET at ~99k cells / 30 cancer types (this session)
- **Criterion 4 (Computational validation)**: MET

## Current Session Work (2026-05-12)

### Kang 2024 Fig 4 — COMPLETE
All 5 NMF compartments integrated with Harmony:
| Compartment | n_obs | n_batches | t_total | OT detection |
|-------------|------:|----------:|--------:|-------------|
| B cells | 9,334 | 662 | 22.4s | — |
| T/NK | 15,689 | 953 | 47.9s | — |
| Mesenchymal | 22,310 | 968 | 76.8s | — |
| Epithelial | 29,084 | 996 | 103.5s | — |
| Myeloid | 22,472 | 1,043 | 356.6s | **p=0.005** ✓ |

**KEY RESULT**: LAMP3+ mregDC OT p=0.005, τ=93.82, Δ=3.337σ, n·π²·Δ²=25.08 (above Theorem 1 floor)

### REAL Motif Parquet — COMPLETE
- workspace/results/kang2024_fig4/myl_real_motifs.parquet
- M-KANG-MYL-001: LAMP3+ mregDC, n=225, π=0.01, Δ=3.337σ, OT p=0.005
- Tumor Biology + Immunology 24h annotation SLA started

### Melanoma Slice — IN PROGRESS
- Cutaneous (mel_GSE200218): 28,291 cells, 5 batches, OT detection running (CPU)
- Uveal cohorts: queued as secondary run

### Git Status
- agent/computational_biology@e2bd7a1 (clean, CB files only)
- .gitignore added to prevent stray dirs

## Pending Work (this cycle)
1. Cutaneous melanoma OT result → notify Immunology
2. Per-cancer-type D-column for myl motif (Tumor Biology request)
3. Scanorama+scVI on Kang myl (scale-dependent inversion confirmation for ML)
4. Baron/Muraro/Segerstolpe/Wang pancreas download (§08 un-caveat)
5. Caushi 2021 NSCLC download (GSE173351)

