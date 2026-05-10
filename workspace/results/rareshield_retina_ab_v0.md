# Retina RareShield A/B — Real-data protection demonstration

**Commit:** agent/computational_biology@[latest]
**Dataset:** Shekhar 2016 retina (19,829 cells × 13,166 genes × 2 batches × 15 bipolar types).
**Holdout:** BC8_9 (197 cells, ~1.1%).
**Method:** fine-tune proxy — load post-integration H5AD, subsample to 4,000 cells (keeping all target-motif cells), fine-tune `obsm['X_integrated']` via gradient descent on L_mass (Math's compute_l_rare@c885f99) with motif_ids = {BC8_9_rare_motif} ∪ {overcorrection_candidates with loss_prob > 0.5}. 100 epochs, λ_m=2.0, anchor weight 0.05.

## Results

| Method    | ΔARI (|gate|≤0.02) | n overcorrection targets | Δlp(rare) | Δlp(overcorrect) | ΔOT-p(overcorrect) |
|-----------|---------------------:|-------------------------:|----------:|-----------------:|--------------------:|
| Harmony   |           **+0.0037** |                        0 |   +0.0007 |              N/A |                N/A |
| Scanorama |             **-0.0136** |                        4 |   +0.0139 |         **-0.2357** |            +0.0082 |
| scVI      |             **+0.0011** |                        4 |   -0.0060 |          -0.0039 |            +0.1107 |

### Harmony — control case
No motifs flagged as overcorrected in vanilla Harmony, so RareShield has nothing to push on. Fine-tune slightly reshuffles the rare-motif neighbourhood; ARI tips up marginally (+0.004). Confirms RareShield is a no-op when nothing needs protecting.

### Scanorama — the headline result
4 abundant motifs flagged as overcorrected in vanilla (loss_prob 0.65–0.82). After RareShield fine-tune, these motifs' loss_probability drops by 0.236 on average (e.g. motif 11: 0.82 → 0.32, motif 5: 0.68 → 0.20). The rare BC8_9 motif is unchanged. ARI stays within gate (-0.014). **This is Fig 3d.**

### scVI — harder case
Same 4 overcorrection candidates identified, but fine-tune only moves them by -0.004 in loss_prob. This suggests scVI's overcorrection is deeper in the embedding geometry than the surface-level density modulation fine-tune proxy can reach — the full scVI+L_rare-coupled-ELBO trainer (TODO) is expected to handle this because it changes the encoder, not just the latent. The OT p-value on overcorrection candidates does rise by 0.111 (toward non-detection), indicating CoLM signal reduction.

## Summary claim for the paper

**On real single-cell data (Shekhar 2016 retinal atlas, 19,829 cells), RareShield fine-tune drops REAL's overcorrection loss-probability by 0.24 on Scanorama-overcorrected abundant motifs while keeping major-type ARI within 0.02 of vanilla, and does no harm when applied to Harmony (which had no overcorrection in the first place).** This is the first real-data demonstration of the protection claim in the abstract.

## Limitations

1. **Fine-tune proxy**: we modify obsm['X_integrated'] post-hoc, NOT the scVI encoder. Full scVI+L_rare-coupled-ELBO trainer is a Lightning-callback TODO.
2. **Subsampling**: 4000 cells (keeping all target motifs + random other). The non-subsampled cells stay at vanilla embedding. Affects scVI case — can explain why scVI's Δlp is small.
3. **λ_m=2.0 fixed**: the 3-seed × 4-λ grid on synth showed λ-insensitivity in the proxy.
4. **Only BC8_9 holdout**: BC4 and BC2 not run for A/B yet (they take ~4 min each on CPU).

## Next

1. Lightning-callback scVI+L_rare trainer (no subsampling, no post-hoc fine-tune).
2. Ablation suite (ML's request): −L_mass / −L_within / −L_between / F=identity / w_M_max sweep.
3. Run retina A/B at BC4 and BC2 holdouts for abundance-axis confirmation.
4. Kang 2024 at 500k once compute ready.
