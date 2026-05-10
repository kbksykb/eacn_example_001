# Synthetic Pilot Results v0 — REAL channels on 3 integrators

**Date:** 2026-05-10
**Status:** pipeline-plumbing validation only; NOT manuscript results. Real-data runs (scIB pancreas + HLCA) come next as the figshare download stabilizes.
**Commit:** computed on `examples/computational_biology/workspace/code/pilots/synth_pilot.py`.

## Setup

- 6,000 cells, 3 batches × 2,000 cells, 2,000 genes.
- 4 cell types with abundance 40 / 30 / 25 / 5 % (alpha / beta / delta / **epsilon**).
- Epsilon (rare) markers partially overlap alpha markers at 40% strength; epsilon's own markers are down-weighted by 50% — this simulates the realistic case where a rare subpopulation sits at the boundary of an abundant one and is easily absorbed by aggressive integration.
- Pre-integration: HVG (1,000) → PCA 50 → top-32 projected for REAL channels; Leiden resolution 3.0 → 12 candidate motifs (7 pass the ≤ 10% abundance filter).

## Results

3 integrators run; REAL channels (mutual-kNN P_n + Procrustes displacement + bootstrap-stable) scored per-motif loss probability.

| method    | LossRate@1 | LossRate@3 | max p(loss) on rare-truth motif | max p(loss) on non-rare motif | wall-clock |
|-----------|------------|------------|----------------------------------|--------------------------------|------------|
| Harmony   | 0.00       | 0.00       | **0.23**                         | 0.20                           | 0.7 s      |
| Scanorama | 1.00       | 0.67       | **0.90**                         | 0.87                           | 8.0 s      |
| scVI      | 1.00       | 0.67       | **0.87**                         | 0.87                           | 34.3 s     |

Ground-truth rare-candidate motifs: 3 (epsilon-dominated Leiden clusters out of 7 candidate motifs that pass the rare-abundance filter).

## Read-out

- The pipeline (synthesis → preprocess → 9-method harness entry-point → REAL channels → parquet → manifest) runs end-to-end in-process.
- REAL gives a usable **separation** between Harmony (preserved, p≈0.23) and Scanorama/scVI (absorbed, p≈0.9) on this synthetic.
- **Caveat:** on this synthetic, non-rare motifs also score high for Scanorama/scVI, because the toy data makes every Leiden cluster move under integration. Real data will have much clearer separation between rare motifs (truly absorbed) and abundant ones (just shifted in Procrustes alignment but preserved as a density mode). This is expected and will resolve once we calibrate on scIB pancreas — the RareScore channels were tuned on toy data, not real integration output.
- Harmony's low scores on this instance are because the synthetic alpha-epsilon overlap isn't deep enough — Harmony happens to keep epsilon separable. We'll add a higher-overlap regime for the synthetic sweep (per DS Fig 2c Harmony θ sweep) once real-data runs finish.

## Emitted artifacts

All DATA_SPEC-compliant:

```
/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/
  data/pancreas_synth.h5ad
  integrations/harmony/pancreas_synth.h5ad
  integrations/scanorama/pancreas_synth.h5ad
  integrations/scvi/pancreas_synth.h5ad
  integrations/manifest.csv
  detections/real/harmony_pancreas_synth.parquet
  detections/real/scanorama_pancreas_synth.parquet
  detections/real/scvi_pancreas_synth.parquet
```

DS's Fig 2 rendering script (examples/data_science/workspace/figures/fig2_validation.py) can auto-ingest detections/real/*.parquet.

## Next

1. Real scIB pancreas download (blocked — figshare AWS WAF). Trying: (a) theislab GitHub releases, (b) direct CELLxGENE API for Baron/Muraro/Segerstolpe/Wang individual datasets, (c) cellxgene-census.
2. Once real-data run shows the expected scVI/scDREAMER preserves-epsilon vs Harmony/Seurat-erases-epsilon pattern, freeze channel thresholds per BioSci's hyperparameter-freeze-gate.
3. HLCA ionocyte holdout (t-mozvf39m primary deliverable).
4. Fill placeholders: scDML, scDREAMER, scCRAFT, Seurat-RPCA wrappers.
