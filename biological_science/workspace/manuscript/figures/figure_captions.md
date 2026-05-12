# Figure captions — main figures (Data Science draft)

These captions follow the manuscript nomenclature proposed in workspace/05_naming_unification.md:
- **REAL** — label-free detection framework (formerly LRS).
- **L(f; w)** — the per-motif survival score (formerly LRS(w)); **L(f)** = mean over the witness set **W**.
- **RareShield** — the motif-preserving integration method (ML agent's proposal).
- **motif** — cross-batch-reproducible density anchor seed (formerly "witness / seed").

Captions below are drafts — Biological Science owns final edits.

---

## Fig 1 | Detecting and protecting previously-unannotated rare subpopulations is a label-free problem.

**a**, Illustration of the core problem. Two batches of scRNA-seq data contain the same three abundant populations (blue circles, batch 1; green squares, batch 2) plus a shared rare subpopulation (orange). A standard integration operator removes the global batch shift but also diffuses mass in the embedding space, absorbing the rare subpopulation into its nearest abundant neighbor. Because this population carries no pre-existing label, it disappears silently from every downstream analysis.

**b**, The **REAL** pipeline for label-free detection. Step 1: compute a per-batch density anchor $R_b(x) = \rho_b(x) / \widetilde{\rho}_b(x)$, where $\rho_b$ is the $k$-NN density of cell $x$ in its own batch and $\widetilde{\rho}_b$ is the smoothed density over a larger $K$-NN neighborhood. Step 2: form a cross-batch motif witness set $W$ by connecting per-batch seed clusters via mutual nearest neighbors on highly-variable-gene (HVG) centroid cosine. Step 3: for each motif $w \in W$, compute a survival score $L(f; w) = \mathrm{HM}(P_d, P_t, P_n, 1-B)$ combining density persistence, topological persistence, neighborhood integrity, and compositional bleed. The whole pipeline scales to 4.9 M cells via FAISS-GPU $k$-NN, cugraph connected components, and witness-complex persistent homology.

**c**, The **RareShield** integration objective. A scVI-style encoder $f_\theta$ is trained with the standard reconstruction and batch-adversarial losses on non-motif cells, plus two motif-preserving regularizers: $\mathcal{L}_\mathrm{anchor}$ fixes motif centroids across batches, and $\mathcal{L}_\mathrm{mass}$ enforces mass-preservation on motif cell clouds via an entropy-regularized optimal-transport term. The adversary is blinded to motif cells, so it cannot drive overcorrection at the cost of rare-population geometry.

**d**, Toy demonstration on synthetic two-batch data with one injected rare cluster (20 cells out of 620, 3%). A standard integrator (strength 0.7) yields $L = 0.18$ — the rare cluster is eliminated. A RareShield-style integrator with matched strength yields $L = 0.84$ — the rare cluster is retained while major populations are still aligned.

Schematic; no labels are used by either step of REAL.

---

## Fig 2 | REAL detects loss of held-out rare populations without using their labels.

**a**, Pancreatic islet atlas (Baron, Muraro, Segerstolpe, Wang; $n \approx 25\,000$ cells) pre- and post-integration (Harmony) UMAPs, colored by REAL motif score. Epsilon, Schwann, and activated-stellate cells — which are the classical "test cases" for rare-subpopulation preservation — all appear as motifs in the pre-integration view; three of them are erased post-Harmony and are flagged by REAL.
**b**, Mean-z-score marker heatmap over REAL-flagged motifs versus published epsilon / Schwann / activated-stellate signatures; rows are motifs ordered by $L$ (ascending).
**c**, Dose-response: $L$ as a function of Harmony $\theta$, showing a monotone dependence on regularization strength.
**d**, Same analysis for pulmonary ionocytes (HLCA core, Sikkema 2023; $n \approx 2.4$ M) and macaque retinal OFFx cells (Peng 2019).

Data will be filled in once Comp-Bio delivers `shared/integrations/<method>/<dataset>.h5ad` per the data spec in workspace/manuscript/figures/DATA_SPEC.md.

---

## Fig 3 | RareShield protects motifs without compromising major-type integration.

**a**, Per-method ranking on the REAL benchmark: box-plot of $L$ across motifs, for Harmony, Seurat-RPCA, Scanorama, BBKNN, scVI, scANVI, scDML, scDREAMER, scCRAFT, and RareShield; datasets pooled.
**b**, Classic scIB / scIB-E metrics on the same integrations: ARI, NMI, ASW, graph connectivity, iLISI, cLISI. RareShield matches baselines on major-type metrics while dominating on rare-type metrics.
**c**, Synthetic injection: $L$ vs. rare-population frequency (0.05 %–2 %) and vs. transcriptomic separation $\Delta\mu$; gives the detection-protection envelope.
**d**, Ablation: remove $\mathcal{L}_\mathrm{anchor}$, remove $\mathcal{L}_\mathrm{mass}$, replace admissible field $F$ with identity, vary $k$.
**e**, Runtime / peak memory on synthetic 100 k / 1 M / 4.9 M-cell datasets.
**f**, Per-method UMAP side-by-side on pancreas, highlighting preserved vs. absorbed rare clusters.

---

## Fig 4 | Pan-cancer immune atlas — REAL surfaces previously unannotated motifs at scale.

Kang et al. 2024 pan-cancer immune atlas, 4.9 M cells, 30 cancer types.
**a**, Global UMAP after standard integration (Harmony / scVI), colored by cancer type (hex-binned for legibility at scale).
**b**, Global UMAP after RareShield-integrated embedding, with REAL-flagged motifs highlighted.
**c**, Volcano: per-motif reproducibility (fraction of cancer types in which the motif is witnessed) vs. per-motif $L$ under Harmony; top-left quadrant = high-reproducibility motifs lost by Harmony.
**d**, Top-10 motif signature heatmap, annotated by Tumor Biology / Immunology against the rare-subpopulation registry (TCF7+CD8, CD103+CD39+ TRM, CCR8+ eTreg, LAMP3+ mregDC, TREM2+/FOLR2+/SPP1+ macs, CAF sub-phenotypes, EMT-hybrid, drug-tolerant persisters, Tip EC).
**e**, Cross-atlas replication matrix: motif reproducibility between Kang 2024 and a second atlas (selected by Tumor Biology, candidates: Gao 2024 / Zhang 2020 / Bi 2021).

---

## Fig 5 | Case study: a novel motif surfaces in one cancer type.

One cancer slice (selected by Tumor Biology based on §8 priority ranking — melanoma first). Panels: (a) UMAP with the new motif highlighted, (b) marker heatmap + biological annotation by Immunology, (c) cohort and batch-of-origin consistency, (d) differential expression vs. nearest canonical subset, (e) discussion hooks (therapy implications, cohort associations).

---

## Extended Data
- ED1 REAL seed-enumeration robustness: sensitivity to $k$, $K$, $\tau$, $\sigma_\mathrm{thresh}$, minimum witness count.
- ED2 Failure modes where REAL mis-detects (doublet clusters, ambient-RNA artifacts, dissociation-stressed clusters).
- ED3 All 10 integrators, full panel of metrics, all datasets.
- ED4 Pan-cancer runtime / memory profiling.
- ED5 Motif stability across random seeds, subsampling, and batch-balance.
- ED6 Cross-witness ≥2 vs ≥1 comparison.
- ED7 RareShield ablations — per-loss, per-backbone (scVI / scANVI / Harmony-OT).
- ED8 Comparison to CellANOVA-style post-hoc signal recovery.
- ED9 Evaluation on cross-species retina (macaque / human).
- ED10 Sensitivity to clustering resolution at the motif-enumeration step.

---

Captions will be trimmed for Nature length conventions once we have numerical results. Everything above is a working draft for Biological Science to reshape.
