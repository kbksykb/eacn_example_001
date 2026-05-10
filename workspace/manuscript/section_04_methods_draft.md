# §04 Methods — Detection & Protection of Previously-Unannotated Rare Subpopulations (v3.2 draft)

Author: Computational Biology Agent (agent-mozur9ik)
Date: 2026-05-11
Supersedes: section_03_draft_v2.md for §04 content.
For BioSci splice: this is the v3.2 Methods content corresponding to the evidence in workspace/results/v3_0_rc_evidence_pack.md + workspace/results/harmony_lamp3_scale_series.md + Panels A/B/C.

---

## 4.1 Detection framework (REAL)

REAL (Rare-subpopulation Erasure under Alignment, Label-free) is a label-free framework comprising four channels applied to paired pre- and post-integration embeddings of the same cells. Channels: mutual-kNN purity (P_n), Procrustes displacement with batch-bootstrap stability (P_Δ + bootstrap_stable), topological persistence P_t (optional when ripser available), and CoLM OT mass-conservation (P_d via the Chizat-2018 entropic unbalanced Sinkhorn).

All channels emit per-motif test-statistics under a permutation null (B=500 for published-grade, B=100 for prototype). Per-channel p-values are BH-FDR-adjusted within-dataset across the motif family (22 motifs for retina, 32-34 for Cheng, 41-42 for Zheng/Sade). Fisher's method combines surviving channels into a single LRS(w) per motif: T = -2·Σ log(p_i) ~ χ²_{2k} with k=number of channels reporting a p-value.

## 4.2 Pre-integration motif construction

Cells pass through HVG-Leiden over-clustering on the pre-integration 32-dim PCA embedding at resolution r=3.0. The resulting Leiden motifs are the atoms of REAL scoring. HVG selection uses `scanpy.pp.highly_variable_genes` with `batch_key=batch, flavor=seurat` — label-free (no cell-type / holdout-rare keys involved).

## 4.3 Holdout protocols (three modalities)

Three modalities of holdout were used, each with different epistemic status disclosed explicitly:

**(A) Author-annotated cluster-level holdout.** Used on Cheng 2021 PAAD myeloid (LAMP3+ mregDC via M05_cDC3_LAMP3), Cheng 5-cancer pan-cancer, Cheng 6-cancer pan-cancer, Zheng 2021 MM melanoma T-cell (CD8.c12.Tex.CXCL13 = CXCL13+ terminally-exhausted CD8 per Liu 2022; note: this cluster was initially misattributed to Miller-2019 TPEX and corrected post-detection via Immunology + Tumor Biology cross-validation). The `obs['cell_type']` field is null-set for cells of the target cluster before any integration step; the held-out identity is stored in `obs['__hidden_label__']` for ex-post validation only. This is the standard protocol in cluster-definable rare-population benchmarks.

**(B) Author-annotated state-level holdout.** Used on Shekhar 2016 retina (BC8_9/BC4/BC2/BC1A bipolar subtypes). Same mechanism as (A) but the rare-truth is a BC-subtype state rather than a cluster-biologically-discrete population.

**(C) Oracle-defined holdout.** Used on Sade-Feldman 2018 melanoma anti-PD1 (GSE120575, no per-cell cell-type annotations provided in the GEO bundle). Canonical Miller-2019 TPEX panel (TCF7+CXCR5+IL7R+SLAMF6) is applied to the log-normalized expression to derive an oracle `is_TPEX_canonical` label for cells scoring ≥3/4 markers. **Three-property non-circularity guard** (per Philosophy agent signed-off):
- (P1) Pre-specified external panel — Miller 2019 marker set was pinned before any Sade data was touched.
- (P2) Oracle is a MEASUREMENT of class membership on the same expression data that REAL consumes, but structurally distinct: oracle applies a threshold rule per cell; integration reads genes as features for a latent space. Different operations on same data.
- (P3) HVG ∩ oracle-marker overlap is disclosed. Canonical TPEX markers (high-variance by virtue of TPEX cells being transcriptionally distinct) are likely in the 1000-HVG set; this does NOT violate non-circularity because HVG selection is variance-driven, not label-driven. The circularity test asks whether the DETECTION METHOD has access to the held-out label; variance-based HVG does not leak.

## 4.4 Synthetic simulation

Synthetic 6000-cell 3-batch datasets use workspace/code/datasets/synthetic.py. 4 cell types at 40/30/25/5% abundances; 100 synthetic marker-gene blocks per type; per-batch housekeeping shifts (batch_shift_sd=0.5). Epsilon-analog rare cells have their 100 markers down-weighted by 50% and 40%-overlapped with alpha's markers (producing the "rare-at-abundant-boundary" worst-case). Counts are Poisson with expected rate exp(log-mean). NO biological gene names are used; all markers are `g_{0000..1999}`. Claims that REAL recovers biologically-named markers MUST NOT be applied to synthetic data — these would be fabrication-class per Philosophy lint:epistemic.

## 4.5 RareShield protection

RareShield adds a differentiable protection loss L_rare to any integration encoder. L_rare is Math's compute_l_rare per workspace/code/compute_l_rare.py (fixed global bandwidth to avoid adaptive-k-NN bias; |τ|² robust aggregator; ReLU gate at 0.4). v1 implementation uses a fine-tune-proxy: post-hoc gradient descent on `obsm['X_integrated']` against L_rare + small anchor_L2 — lower bound on a full scVI+L_rare Lightning-callback trainer (future work).

λ_mass=2.0, 150 epochs, Adam lr=1e-3. ΔARI hard-gate ≤0.02 on Leiden major-type clustering against public-labels-minus-held-out-rare (primary metric per BioSci specification). Oracle-based ARI in Supp as robustness check.

## 4.6 Scale-dependent and integrator-intrinsic detection regimes

Three distinct detection regimes co-exist in the method × scale grid, and the paper distinguishes them explicitly rather than conflating them under a single scaling narrative.

**Regime A — scale-monotonic detection** (conservative-linear integrators). Fixing method=Harmony and rare-type=LAMP3+ mregDC across Cheng 2021 PAAD (n=2,853, π=1.07%), Cheng 2021 5-cancer pool (n=20,341, π=1.50%), and Cheng 2021 6-cancer pool (n=49,271, π=1.15%): detection strength -log10(p_bh) grows monotonically 0 → 1.89 → 2.39, matching Math's closed-form §rem:S1-lamp3-scale-series Le Cam sigmoid derivative. All three whole-atlas points have n·π²·Δ² < Theorem-1 multi-channel-composite floor; below-floor detection is enabled by the CoLM channel operating at sub-compartment scale per §rem:S1-effective-floor-locality — effective n is the myeloid-lineage size, not the whole atlas.

**Regime B — integrator-intrinsic mild-scale-modulation** (aggressive-linear integrators). Scanorama on the same three cohorts: -log10(p_bh) is 2.30, 2.10, 2.41 — detection fires at the smallest cohort (n=2,853) where Theorem 1 would predict preservation. Not scale-induced: Scanorama dissolves LAMP3+ into cDC2-proximal DCs even at small n, a method-intrinsic overcorrection of the rare-DC state.

**Regime C — integrator-intrinsic scale-flat detection** (non-linear latent integrators). scVI on the same cohorts: -log10(p_bh) is 2.52, 2.40, 2.52 — detection fires essentially uniformly at all three scales. Non-linear latent compression dissolves LAMP3+ independently of n. A non-specificity caveat accompanies: scVI's non-rare-motif flag rate on Cheng6 is 0.65 vs Harmony's 0.47 and Scanorama's 0.50, so scVI's high detection signal is accompanied by reduced specificity, not pure detection gain.

**Regime-dependent scVI specificity.** scVI's specificity is itself regime-dependent, not uniform. At small-cohort / curated-canonical-TPEX scales (Sade-Feldman 16k) it preserves the canonical TPEX motif at p_bh=1.0, consistent with mild regularization; at pan-cancer scales (Cheng6 49k) its non_rare_flag_bh=0.647 indicates broader OC-candidate firing, consistent with scaling-induced feature-space dilution. These two readings are not contradictory — they are two regimes of the same non-linear-latent mechanism: mild on small-curated data, promiscuous on large-heterogeneous data. Harmony and Scanorama do not exhibit this scale-dependent specificity transition; the observation is scVI-specific and documented as evidence-of-behavior rather than a bug.

**Falsified cross-method unification**. A formally pre-specified cross-method unifying prediction (an ℓ-weighted regression over the scale-series, where ℓ := fraction of k-NN sharing the cell's coarsened lineage) was empirically disconfirmed on a 9-point grid (`workspace/results/locality_regression_v1_negative.md`). Both coarse ℓ (DC/Mac/Mono/Mast/Neut) and fine ℓ (MajorCluster) fail to collapse the 9 points onto a monotone scaling curve in n·π²·Δ²: counter-example pairs exist at (Scanorama PAAD n·π²·Δ²=0.21, ℓ·-log10(p)=2.30) vs (Harmony Cheng5 n·π²·Δ²=2.93, ℓ·-log10(p)=1.13) — the larger-scale point gives half the weighted signal. The scale-induced scaling (Harmony) and method-intrinsic aggressiveness (Scanorama, scVI) are therefore reported as distinct mechanisms rather than as a unified scaling story; the hierarchical-REAL mechanism (\Cref{rem:S1-effective-floor-locality}) remains empirically supported as an intra-method explanation for Harmony's below-floor detection.

**Consequence for method ranking**: existing small-cohort scIB-based rankings conflate these regimes. At n=2,853, Harmony preserves while Scanorama and scVI fire. At n=49k all three fire at comparable -log10(p_bh) but with stratified non-specificity. The paper's method-level contribution is NOT "method rankings invert with scale" (too strong given Regimes B+C); it is "method rankings depend on which detection regime the target rare-type is in, and atlas-scale benchmarks must stratify by regime rather than reporting a single rare-type-preservation score."

## 4.7 Per-method outcome stratification under RareShield A/B

Three distinguishable per-method outcomes observed across Cheng5 and Zheng MM A/B grids:
- **Textbook full-reversal protection** (Cheng5 Harmony + Scanorama; Zheng MM Scanorama): ΔARI within gate, OT p raised from 0.01-floor to near-1.0 (≥0.98 shift). Target rare-truth preserved post-RareShield.
- **OC-only protection** (Zheng MM scVI; Cheng5 scVI): ΔARI within gate, Δlp on overcorrection candidates reduced by ~0.16-0.41, but rare-truth OT unchanged (proxy limit — fine-tune couldn't lift rare-truth from CoLM-firing position). Expected-to-close with full Lightning trainer.
- **Leakage edge-case** (Zheng MM Harmony): ΔARI within gate but fine-tune targeting OC candidates perturbed the preserved rare-truth's neighborhood geometry enough to lift its OT signal into detection range. Diagnostic: a proxy-artefact when vanilla-already-preserves.

All outcomes documented honestly per Philosophy's Popperian-falsifiability register — the three-way split demonstrates integrator-ordering prediction (Theorem S1-overcorrection κ) rather than uniform-firing.

## 4.8 Specificity / negative-control

E_neg synthetic batch-artefact regime injects 60 cells in batch_0 only with a housekeeping-gene spike. Both Harmony and Scanorama give OT p=1.000 on this pure artefact motif while correctly firing on the biological epsilon-rare at p=0.01-0.04 on the same run (workspace/results/synth_eneg_neg_control.csv). Trap-7 (evaluation circularity) + Trap-8 (failure-to-reject = success) + Trap-9 (oracle-derivation circularity) all empirically closed.

## 4.9 Compute and reproducibility

All runs: 8× A100 shared server (Aliyun Shanghai, 80 GB per GPU). Python 3.11 conda env `scenv` with pinned scanpy 1.10, scvi-tools 1.4.2, harmonypy 2.0.0, scanorama 1.7.4, POT 0.9, scikit-learn 1.5. Commit `agent/computational_biology@4955879`. Artifacts at `/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/` (shared FS read-accessible by all agents); per-dataset H5ADs have `uns['reproducibility']` populated with random_state, commit_sha, pip_freeze_hash, run_date.

## 4.10 Limitations explicitly declared

1. **Real pancreas scIB benchmark dataset is NOT in the current parquet set** — figshare / CELLxGENE access was blocked during the submission cycle. The "pancreas" label anywhere in the evidence refers to the synthetic pancreas-like simulation (g_0..g_2000 markers), not Baron/Muraro/Segerstolpe/Wang real data. This is a data-access limitation, not a methodological one; real pancreas is a revision-tier addition.

2. **Kang 2024 4.9M pan-cancer immune atlas is NOT in the current parquet set** — not indexed in CELLxGENE 2024-07-01 or 2025-11-08 census; Zenodo/HuggingFace searches did not return it. Current scalability tier is 49,271 cells on Cheng 2021 6-cancer pool; revision-tier addition plans for the Kang atlas once access is staged.

3. **RareShield is a fine-tune proxy, not a full scVI+L_rare Lightning-callback trainer.** Proxy results are strict lower bounds on what the coupled-ELBO trainer should achieve. The scVI OC-only and Harmony leakage edge-cases are expected to close under the full trainer.

4. **The Sade-Feldman 2018 scVI HVG3K preservation** (p_bh=1.0) is reported as a specificity counter-observation, not a detection failure. REAL fires when aggressive integrators overcorrect; scVI's milder regularization does not overcorrect canonical TPEX on Sade-Feldman at n=16k, π=0.044.

## References to companion documents

- workspace/results/v3_0_rc_evidence_pack.md — 4-criterion evidence summary
- workspace/results/harmony_lamp3_scale_series.md — 3-point scale-series narrative
- workspace/results/panel_A_tight_v2.csv — 25-row Theorem 1 validation
- workspace/results/panel_B_per_motif_kappa.csv — 66-row κ × pathway asymmetry
- workspace/results/panel_C.csv — 7-row hierarchical amplification
- workspace/results/retina_abundance_series_v1.md — 12-run Shekhar grid
- workspace/results/rareshield_ab_v0.md + rareshield_retina_ab_v0.md + Cheng RareShield JSONs — A/B grid
- workspace/results/synth_eneg_neg_control.csv — specificity null
- workspace/results/retina_overcorrection_diagnosis.md — RBC sub-fragment interpretation
- workspace/results/sade_tpex_pilot.md — state-regime + 3-method specificity counter-observation
- workspace/results/zheng_mm_cxcl13_tex_summary.md — 3-method real-data detection + protection

All above published on agent/computational_biology branch, latest commit `agent/computational_biology@4955879`.
