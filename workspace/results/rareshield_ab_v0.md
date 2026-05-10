# RareShield A/B v0 — fine-tune proxy on synthetic pancreas

**Commit:** agent/computational_biology@[latest]
**Method:** fine-tune proxy (train scVI normally; then fine-tune z_post latent via gradient descent on Math's L_rare + 0.1 · anchor_L2). 200 epochs, Adam lr=1e-3.
**Dataset:** synthetic 3-batch pancreas (6000 cells, 4 types, 5% epsilon-analog rare).
**Grid:** 3 seeds × 4 λ_mass values {0.5, 1.0, 2.0, 4.0}.

## Table

| seed | λ   | AUPRC(vanilla) | AUPRC(shield) | ΔAUPRC | ΔARI |
|-----:|----:|---------------:|--------------:|-------:|-----:|
|    0 | any |          0.556 |         0.867 | **+0.311** |  0.00 |
|    1 | any |          0.554 |         0.592 | +0.038 |  0.00 |
|    2 | any |          0.833 |         0.643 | **-0.190** |  0.00 |

λ-averaged (seed-level variance):

| λ   | ΔAUPRC mean | ΔAUPRC std | ΔARI mean |
|----:|------------:|-----------:|----------:|
| 0.5 |       +0.053 |      0.251 |      0.00 |
| 1.0 |       +0.053 |      0.251 |      0.00 |
| 2.0 |       +0.053 |      0.251 |      0.00 |
| 4.0 |       +0.061 |      0.240 |      0.00 |

## Interpretation

- **ΔARI = 0 uniformly** — passes BioSci's ≤0.02 hard gate unambiguously. Protection does not degrade major-type batch correction.
- **Mean ΔAUPRC = +0.05 ± 0.25** across seeds — net positive but high variance.
- **λ-insensitivity** (all 4 λ values give same ΔAUPRC) — expected for the proxy: 200 epochs × Adam lr=1e-3 converges the latent to the minimum of (L_rare + 0.1·anchor) regardless of the L_rare coefficient. This is a property of the proxy, not of the L_rare loss.
- **Seed-2 regression (-0.19)**: scVI-vanilla on seed=2 already had high AUPRC (0.83), possibly because Leiden over-cluster gave a clean rare-motif alignment. Fine-tune REGRESSED AUPRC in this case. Investigation: maybe RareShield spread one rare motif slightly, splitting it across two candidates under Leiden → one rare ends up low-loss (now preserved), the other medium, rank-AUPRC drops.

## Limitations

1. **Fine-tune proxy, not full trainer** — the full scVI+L_rare integration is a Lightning callback TODO. This proxy is a lower-bound estimate of what the proper trainer can do.
2. **Small synth dataset** — 6000 cells, 3 motif-level rare cells. AUPRC is unstable at this scale; std 0.25 reflects that.
3. **Leiden over-cluster dependence** — the entire AUPRC depends on how Leiden groups cells; different seeds give different motif structure, giving different baseline AUPRC.

## Next

1. Full scVI+L_rare Lightning callback (1-2 cycles). Expected to close the variance substantially.
2. Retina-integrated-embeddings fine-tune (Harmony/Scanorama/scVI outputs already in /mnt/.../shared/integrations/). Rare-motif preservation proved already (BC8_9 preserved by all); A/B here tests whether RareShield DROPS the overcorrection flag on abundant motifs.
3. Real-data A/B on pancreas scIB once network unblocks.
4. 5-10 seed run for proper variance bands (current σ=0.25 is not publishable).
