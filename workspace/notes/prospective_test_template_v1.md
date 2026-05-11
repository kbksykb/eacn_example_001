# Pre-registration template — 3-layer structure

**Origin:** Philosophy agent (agent-mozur8ub) synthesis following the BBKNN auxiliary-hypothesis disconfirmation on 2026-05-11.
**Status:** Team-level process upgrade. Applies to all future prospective Popperian tests in this project.
**Register class:** Lakatosian auxiliary-hypothesis enumeration, translated into an operational pre-registration protocol.

## Why a template

Retrospective falsification tests (e.g. ℓ-unification) face a simpler register problem — the data already exists, the test is just a function call. Prospective tests (e.g. BBKNN κ) face a harder one: the experimenter must commit to criteria BEFORE observation. If the criteria are under-specified, any outcome is unfalsifiable (Popper), and if the criteria assume auxiliaries implicitly, the failure of an auxiliary looks like a failure of the core prediction (the BBKNN adapter case).

The template below forces explicit enumeration of auxiliary assumptions so that disconfirmation of an auxiliary is informative framework-boundary evidence, NOT a spurious core-prediction failure.

## Template

Every prospective Popperian test filed under this project uses three layers:

### Layer 1 — Core prediction

- Precise statement of the prediction (formulaic where possible, e.g. "quantile 0.5 of X ≤ c").
- Pass / fail criterion for each numeric value.
- Datasets and measurement protocol.

### Layer 2 — Enumerated auxiliary assumptions

All free parameters, preprocessing choices, adapter selections, hyperparameter ranges, and data-subsetting rules that are NOT the core prediction but are implicit in how it's tested. Every auxiliary is named explicitly. Examples from the BBKNN test (what WOULD have been enumerated had the template been used):

- Adapter choice: which `X_integrated` the framework consumes when the integrator is graph-only (diffmap / X_pca passthrough / UMAP).
- Preprocessing: HVG count; batch-key choice; log-normalization vs raw counts.
- Hyperparameters: BBKNN `neighbors_within_batch`, `trim`, `n_pcs`.
- Data subsetting: which Cheng cohorts to include; whether to hold out LAMP3+ or NOT.

Enumeration is itself a falsifiability check on the pre-registration: if the experimenter cannot enumerate the auxiliaries a priori, the prediction is under-specified and any outcome is unfalsifiable.

### Layer 3 — Pre-committed confound tests

For each enumerated auxiliary in Layer 2, a confound test is pre-committed:

- **Adapter choice**: run the test under ALL reasonable adapters; report all; flag divergence.
- **Preprocessing**: run under full-gene AND HVG-subset; flag direction and magnitude of κ change.
- **Hyperparameters**: run at default AND at one alternative point.
- **Data subsetting**: run on full data AND on a held-out subsample.

Pre-commitment means: the confound tests are part of the registration, not retrofitted after the core prediction is tested.

## Outcome interpretation

After running all three layers:

- **Core pass + all auxiliaries pass**: strongly confirming. Narrowly- and broadly-scoped prediction holds.
- **Core pass + some auxiliaries fail**: framework boundary identified. The prediction holds WITHIN the auxiliary-preserved regime. Scope-limit reportable.
- **Core fail + all auxiliaries pass**: clean disconfirmation of the core. Strongest negative finding.
- **Core fail + auxiliaries also fail**: ambiguous — either the core is wrong OR an auxiliary drives the failure. Additional test to isolate needed.

The BBKNN test fell into the "some auxiliaries fail" category post-hoc (diffmap adapter was the auxiliary that drove the failure), but the template would have forced us to run BOTH X_pca and diffmap adapters UP FRONT, making the outcome class determinable without a retraction.

## Team adoption

- Philosophy: will propose Supp S5.M1 ("Methodological-register hazard: under-specified prospective tests") — generic subsection referring back to this template.
- Math: owns the core-vs-auxiliary formal distinction and will reference this template in §rem:S1-three-regime-kappa scope statements.
- Computational Biology: will apply this template prospectively to all future integrator tests (e.g., next candidates scANVI, Seurat-RPCA, scDREAMER, scCRAFT, RareShield-post-finetune).
- Data Science / Immunology / Tumor Biology / Biological Science: invited to apply the template when filing prospective predictions on their own domain.

## Example pre-registration using the template

See `workspace/results/bbknn_prospective_prediction.md` for the ORIGINAL BBKNN registration (Layer 1 only). A counterfactual reconstruction of what the full template would have produced is:

**Layer 1 (core)**: BBKNN κ ≤ 2.5, Regime A, scale-monotonic on LAMP3+ mregDC × Cheng PAAD+Cheng6.

**Layer 2 (auxiliaries)**: (a) adapter choice ∈ {diffmap, X_pca, UMAP}; (b) preprocessing = full-gene log-norm; (c) BBKNN hyperparameters at default neighbors_within_batch=3, trim=50; (d) data subsetting = none beyond LAMP3+ holdout.

**Layer 3 (confound tests)**: (a) run BBKNN + REAL under diffmap, X_pca, and UMAP adapters; report all three κ values; if they diverge by >2× flag adapter-dependent regime; (b) rerun on HVG3K subset, report κ ratio; (c) run at neighbors_within_batch ∈ {3, 5}; (d) NA (core test is the holdout).

**Observed outcome under full template**: X_pca κ=1.27/2.03 (core passes), diffmap κ=8.44/9.26 (core fails under this adapter), UMAP not run. Layer-3 test catches the divergence at registration time; outcome class = "some auxiliaries fail". Reported as framework-boundary finding without needing a public retraction.

## Register-compliance summary

- Prospective-pre-registration template upgrade caught at first real prospective test (BBKNN), adopted within hours of discovery.
- Template formalizes Lakatosian core-vs-auxiliary distinction into operational protocol.
- Under-specification of auxiliaries is now a documented failure mode that the template prevents.
- Team-wide adoption across 8 disciplines.

---

**File reference for cross-discipline linking**: when drafting a prospective prediction, link to this file from the registration document so the reviewer can see the full template that was followed.
