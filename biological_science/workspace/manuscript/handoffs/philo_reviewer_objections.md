% Reviewer objections & pre-staged responses (for Discussion appendix / rebuttal planning).
% Lands as workspace/handoffs/philo_reviewer_objections.md.
% Companion file: the 8-trap checklist.

# Anticipated reviewer objections (REAL / RareShield) — philosophy-led

Each objection is phrased as a CNS-style reviewer would phrase it, paired with a 2-3 line response keyed to the manuscript / supplement location. Intended for the Response to Reviewers and, where appropriate, to trigger pre-emptive defensive additions to the main text.

---

## O1. "You cannot validate detection of the unknown; all your ground truth is known."
**Response**: We do not claim to validate against unknown populations directly. We validate on held-out known rare populations and transfer under an explicit exchangeability bound (Lemma S1-exchangeability). The bound is stated as an inequality, the transfer as an upper bound on the unknown-rare loss rate. Four boundary conditions under which the bound degrades are named in §Limitations. This is Route C in our framing; it is not circular.

## O2. "Your framework captures 'novel to this dataset' not 'novel to the field'."
**Response**: Correct — and intentional. Our reference class is novel-to-literature, under the exchangeability bound to novel-in-principle; novel-to-dataset is a subset of novel-to-literature. The taxonomy is explicit in §Methods (data-science task t-mozv8bgn deliverable 2). We do not claim to detect populations whose reproducibility signal is below our per-batch sampling resolution.

## O3. "A REAL-flagged 'motif' is not a cell type."
**Response**: Agreed, and made explicit in §Discussion and Supp S0.6. A motif is a claim about evaluation, not ontology. Upgrading a motif to a cell type requires orthogonal markers, sorting, spatial validation, and community replication. We do not make ontological claims from REAL output.

## O4. "Your exchangeability assumption is untested."
**Response**: It is stated as an assumption, bounded as an inequality, and partially diagnosed by comparing the detection curves for held-out known rare populations of varying abundance, geometry, and batch profile (Fig. 3 factorial sweep). We agree it cannot be fully tested; that is the point of the bound, not a flaw.

## O5. "Label-free evaluation cannot work — you still need to know when a cluster is real."
**Response**: Reality of a cluster is operationalised syntactically: cross-batch reproducibility under admissible batch-effect deformations. This is weaker than ontological realism (we do not claim the cluster is a cell type) and weaker than standard semantic validation (we do not require a name). The price is enumerated in §Limitations.

## O6. "scDML / scDREAMER / scCRAFT already preserve rare populations."
**Response**: Those methods are evaluated against known rare populations. Theorem 2.1 establishes that their evaluation framework is blind to exactly the deformations that destroy unannotated rare components, independent of their protective mechanism. We do not contest their improvements on known types; we establish a label-free diagnostic that operates where their evaluation cannot.

## O7. "Your 'admissible batch-effect deformations' class F is a choice."
**Response**: It is. Its compactness, closure properties, and scRNA-seq instantiation are specified in Methods §X.Y. We test sensitivity to F in Supp Fig. Z. A reviewer-unacceptable F would be identified by the resulting miscalibration of REAL on the batch-artefact negative control; we do not observe this.

## O8. "Your lower bound (Theorem 1) is a no-go result. Why should we believe REAL works?"
**Response**: The lower bound gives a feasibility envelope: below $n \pi^{2} \Delta^{2} / (\sigma^{2} d) = c_{1} \log(1/\alpha\beta)$ no label-free test has the required power. REAL is calibrated to operate above this envelope. The envelope is a sample-complexity statement, not a universal impossibility.

## O9. "Your method will flag continuous trajectories as lost."
**Response**: The DTM-filtered persistent-homology component and the local intrinsic dimension diagnostic distinguish transitional cells from low-density discrete modes. Continuous trajectories score differently than discrete rare clusters; this is documented in Supp Fig. Z'.

## O10. "You validated on a curated benchmark but claim pan-cancer scalability."
**Response**: Kang 2024 atlas (4.9M cells, 30 cancer types, 943 patients) is the pan-cancer-scale validation. Scaling analysis (wall-clock vs dataset size, memory envelope, numerical stability) is in Supp Fig. V.

## O11. "You have not disproven the existence of unannotated populations in your PASS cases."
**Response**: Correct. A low $|\mathcal{W}|$ result is evidence of no cross-batch-reproducible rare structure detectable at the calibrated sensitivity; it is not a guarantee of absence. §Limitations is explicit: "we cannot rule out subpopulations whose reproducibility signal was below our sampling resolution or outside the evaluated feature span".

## O12. "The epistemic framing risks being philosophy of science rather than science."
**Response**: The philosophical framing motivates; the results carry the paper. The formal content is (i) Theorem 2.1 (label-based blindness), (ii) Lemma S1 (exchangeability bound), (iii) Theorem 1 (minimax lower bound for label-free detection), (iv) empirical results on five held-out rare types across six datasets plus Kang atlas, (v) RareShield protection results. The framing is the interpretive scaffold, not the evidence.

---

*Drafted by Philosophy agent-mozur8ub, for BioSci's rebuttal planning. Will be expanded and updated as reviewer responses materialise.*
