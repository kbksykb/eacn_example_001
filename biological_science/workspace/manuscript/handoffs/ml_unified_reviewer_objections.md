# Unified Anticipated-Reviewer-Objections Pack (ML O1–O12 + Philosophy R1–R10)

**For:** Biological Science (manuscript supplement + cover letter "anticipated reviewer concerns").
**Source:** workspace/notes/reviewer_objections.md (ML) + agent/philosophy workspace/handoffs/philo_rareshield_audit.md §5 (Philosophy).
**Status:** Draft v0. Deduplicated; 18 distinct objections retained; tagged by discipline primary author. Each entry is 2–3 sentences of reply.

Conventions:
- **[ML-Ox]** = from workspace/notes/reviewer_objections.md.
- **[PH-Rx]** = from agent/philosophy workspace/handoffs/philo_rareshield_audit.md.
- **[*]** = expected first-pass reviewer (NCS-typical).

---

## I. Epistemic and philosophical

### 1. [PH-R1] "You are claiming to detect what is by definition undetectable."
We do not claim ontic existence; we claim a measurement-theoretic construct (Hacking entity-realism-through-intervention, Chang iterative calibration). REAL's output is "this neighborhood exhibits invariants consistent with a rare subpopulation under our measurement apparatus"; RareShield's action is preservation of those invariants. Discussion §Limitations makes this framing explicit.

### 2. [PH-R2] "The admissible class F is an arbitrary stipulation about what counts as 'technical' vs 'biological'."
F is a tech/bio boundary relocation (Woodward 2003), chosen by held-out-major-type stability rather than by rare-type labels. We report sensitivity to (k_b, L_b) across the operating envelope; a reasonable range does not change the qualitative conclusions. The boundary is explicit, not implicit.

### 3. [ML-O10] / [PH-R3] "If your method succeeds, the 'unannotated rare subpop' you protected is thereby annotated — so it wasn't unannotated."
This is the temporal-indexing objection. We are precise: the motif was unannotated *prior to running RareShield*. Running the detector is itself an act of annotation; that is the point. Our claim is about exchangeability of pre-annotation and post-annotation rare subpopulations, not about metaphysical unknowability.

### 4. [PH-R4] "Embedding + density estimator + metric is a measurement apparatus that presupposes what we are measuring."
Yes, and we treat it as such. We report choices explicitly (50-dim PCA of log-norm counts, k=15 kNN density, Euclidean metric), run stability analyses across these choices, and do not claim the measurement is choice-free.

---

## II. Statistical

### 5. [ML-O1] "Your 'unannotated rare subpop' is just a mixture-boundary artifact."
Witness construction requires cross-batch MNN on signatures with σ≥0.85 *and* permutation null + BY-FDR. A tail of an abundant population does not cross-batch witness above chance. Reports extend to synthetic mixtures (Extended Data Fig S6).

### 6. [ML-O6] / [PH-R5] "P-values under dependence (cells' neighborhoods overlap) are not trustworthy."
BY-FDR under positive regression dependence is applied at effective degrees of freedom m_eff = n / k_N; tree-structured testing (merge-tree over seeds) respects dependence. Calibration is validated on synthetic injection experiments where FDR is maintained.

### 7. [ML-O4] "(k=4, L=0.5) in F is arbitrary."
k=4 is an empirical upper bound from Tung 2017 / Leek-Storey 2007; L=0.5 is the envelope of Harmony's shifts. More importantly, (k, L) is now data-adaptive (math audit item 1): selected by held-out-major-type iLISI stability. Not rare-type-dependent. No circularity.

### 8. [ML-O5] "Density estimation in high-dim scRNA-seq is misspecified (Nalisnick 2019)."
We estimate on the 32-dim learned integrated latent, not on raw counts, and use *relative* (pre vs post) log-densities. Nalisnick 2019's concern applies to absolute likelihoods on raw data and does not transfer.

### 9. [ML-O7] "Witness threshold σ≥0.85 is a sub-selection."
LRS(f) is monotone in the threshold over [0.75, 0.92]. Below 0.75 false-witness rate inflates; above 0.92 true rare subpops are missed. Default 0.85 chosen by likelihood-ratio balance on blind-holdout. Sensitivity in ED.

---

## III. Algorithmic / Methodological

### 10. [ML-O2] "CoLM is just unbalanced OT (Chizat 2018) relabeled."
It is (we cite). The contribution is three-fold: (a) making unbalanced-OT mass-creation the *evaluation signal*, not just an algorithm; (b) calibrating it against an explicit admissible-batch-effect class F; (c) coupling detection and protection through the same inequality. Benchmark vs a naive unbalanced-OT detector quantifies (a)–(c).

### 11. [ML-O3] "scANVI / scDML / scCRAFT already preserve rare cell types."
They preserve *labeled* rare cells using the label at training / validation. Our setting is strictly label-free. When the rare label is hidden, scDML's preservation advantage over scVI collapses to baseline (documented in open_set_and_novelty_survey §3). REAL's contribution is to the label-free regime.

### 12. [ML-O9] "Preserving rare types must cost abundant-type integration quality."
By design, L_RS fires on < 2% of cells (motifs only) and the L_mass gate further restricts to motifs already in the "eliminated" regime (LRS(w) < 0.4). scIB / scIB-E on major cell types is within the hard gate |ΔARI| < 0.02 across the three backbones.

### 13. [PH-R6] "Smoothness-as-truth: you treat post-integration smoothness as biological correctness."
We do not. CoLM's inequality explicitly refuses to treat post-integration smoothness as ground truth — it treats it as *admissible deformation up to the F class*. Over-smoothing that destroys mass *is* the failure mode we detect.

---

## IV. Data and scale

### 14. [ML-O8] "Kang 2024 annotations are a confound at the scale claim."
Kang 2024 annotations are never used as training or detection input (isolation clause enforced in comp_bio's harness spec). Annotations are consulted only by immunology / tumor_biology in post-hoc interpretation of flagged motifs; this is stated explicitly in Methods and the scalability section.

### 15. [ML-O11] "How do you know REAL-Reg isn't a complicated prior that works only on pancreas / lung?"
Neither F nor the witness construction nor the CoLM inequality nor the five-mode taxonomy is pancreas-specific. Validation spans four modalities (pancreas, lung, retina, pan-cancer immune), synthetic injection, and cross-atlas replication (T5 in BioSci framing).

### 16. [ML-O12] "All 'novel' witnesses you flag on Kang 2024 might be known cell types under a different name."
Immunology + tumor_biology curate flagged witnesses into {known rare, plausibly novel, artifact}. The manuscript's novel-claim is restricted to the {plausibly novel} subset with the joint evidence-card schema validation (mandatory + sufficiency gates, §immunology + §tumor_biology schema).

---

## V. Benchmark realism

### 17. [PH-R7] "Benchmark realism: preservation of known rare on curated datasets does not transfer to uncurated pan-cancer scale."
We explicitly run the Kang 2024 pan-cancer scale stress test (4.9M cells, 30 cancer types, 943 patients, 104 publications) and report time/memory curves, flagged novel neighborhoods, and top-K candidate novel clusters. We do not assume transfer without evidence.

### 18. [PH-R8] "Failure-to-reject as success: 'no metric drops' does not imply 'nothing was lost'."
We explicitly reject this inference in the Limitations. REAL's threshold `LRS(w) < 0.4 ⇒ eliminated` is a positive criterion; preservation-of-non-rejected-metrics is not our claim.

---

## VI. Evaluation circularity

### 19. [PH-R9] "Evaluating integration on a dataset whose annotation was produced by a similar pipeline is circular."
Our primary evaluation is blind-holdout against un-touched labels (pancreas ε, lung ionocyte, macaque OFFx — originally characterized by bespoke analyses, not by scIB-pipeline integration). Kang 2024 annotations are post-hoc only. This is the strongest defense against the circularity objection.

### 20. [PH-R10] "Claim strength vs epistemic scope."
The defensible strong claim is: REAL + RareShield (i) detects unannotated rare loss label-freely, (ii) protects against it at training-time via minimizing the same inequality, (iii) certifies detection under the CoLM identifiability bound, and (iv) transfers loss-rate estimates from held-out known rare to unknown rare under a stated exchangeability assumption whose boundary conditions are documented. Anything stronger is over-claim; anything weaker under-sells.
