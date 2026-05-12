# Reconciliation: Root-Cause Analysis (Path 1) vs. Current Research Methodology

**Author**: Biological Science Agent (`agent-mozuoaqx`)
**Date**: 2026-05-11
**Status**: team-discussion input for manuscript v2.7 prep + Tier-B Kang execution

## 1. Why this reconciliation was requested

The human user flagged `examples/root_cause_analysis_path_1.md` (Triple Absence / Self-Closing Loop) as a path of analysis they want us to compare against the route we have already taken. They also flagged `examples/README.md` Criterion 3 (Scalability on Kang 2024 pan-cancer 4.9M cells) as an immediate priority that is still outstanding.

## 2. What the root-cause path says

The path diagnoses the five-year failure to solve the problem as a **triple absence**, interlocked into a self-closing loop:

1. **No ground truth** → can't define "successfully preserved"
2. **No evaluation framework** → no failure signal
3. **No applicable algorithm** → no protection

Every step waits on the previous. The default escape (external ground truth via simulated data or experimental validation) doesn't generalise: simulated data is circular; annotation-based validation is by construction silent on unannotated populations.

The **unlock**: the path observes that normal batch correction (convergence — same-type cells merged across batches) and pathological rare-subpopulation erasure (dispersion — members of a cluster scattered into unrelated large clusters) leave **structurally different patterns** in the data itself. The signal is accessible without external ground truth; the problem is extracting it.

## 3. How our methodology already aligns

Our pipeline reaches exactly the same entry point but through a different operational vocabulary. The mapping is:

| Root-cause path (Path 1) | Our pipeline | Our agent/file |
|---|---|---|
| "No ground truth" | Label-free surrogate ground truth = set $W$ of cross-batch-witnessed reproducible rare motifs | DS: `agent/data_science:workspace/00_design_lrs_framework.md` §2 |
| "No evaluation framework" | REAL detectability framework with per-seed score LRS($w$) = HM($P_d, P_t, P_n, 1-B$) | manuscript §02 + Math Thm 0 (label-blindness) |
| "No applicable algorithm" | RareShield protective integration with L_RS = $\beta_w(L_{\rm within} + L_{\rm between}) + \lambda_m E_w[\tau^2 \cdot \text{ReLU}(0.4 - \text{LRS})]$ | ML: `agent/machine_learning:workspace/proposals/unified_ML_spec_v3.md`; Math Thm 5 identifiability |

The **convergence-vs-dispersion distinction** Path 1 names as the unlock is literally the object our $P_n$ (neighborhood integrity) and $B$ (compositional bleed) components measure:

- Normal convergence: $P_n(w)$ high (cells of a motif keep their pre-integration MNN neighbours); $B(w)$ low (non-motif neighbours not dominated by a single adjacent abundant cluster).
- Pathological dispersion: $P_n(w)$ low; $B(w)$ high (motif cells absorbed by a nearby abundant population).

Math's **Theorem S1-overcorrection** (agent/mathematics) generalises the same observation: REAL's CoLM statistic fires on any motif (rare or abundant) when integration displaces fraction $\kappa_w$ through a low-density stratum — "dispersion" in Path 1's vocabulary.

## 4. What Path 1 adds that our paper does not yet state explicitly

Four items are worth lifting into the manuscript:

1. **The self-closing-loop diagnosis**. Our Introduction argues "label-based evaluation is blind to unannotated loss" (Math Thm 0). Path 1 sharpens this to a three-step loop that explains *why* five years of effort didn't break it: the motivation gap, not just the technical gap. This is a two-sentence add in §01 (between current paragraphs 3 and 4).

2. **Convergence-vs-dispersion as the named mechanism**. Our §02 defines $P_d, P_t, P_n, B$ technically but never names what they collectively measure in plain language. Path 1's framing — "convergence vs. dispersion leaves structurally different patterns" — is exactly the plain-language statement. One sentence in §02 Intuition would close that gap.

3. **Circular-reasoning critique of simulated-data validation**. Path 1 calls out "you know it exists because you built it" about simulated benchmarks. Our Methods already frames synth as sensitivity and real as specificity (per Immunology's ladder), but the word "circular" does not appear in the Limitations. Adding it strengthens Philo's trap #1 / #6 disclosure.

4. **"The signal exists in the data; no one designed the framework to extract it"**. This is an abstract-worthy sentence. A variant is already in our abstract ("We close this gap with REAL...") but Path 1's more direct form — *the signal exists; a framework was missing* — is punchier. Optional one-line polish.

## 5. Where Path 1 and our paper differ

Two differences are substantive and deliberate:

- Path 1 treats the unlock as **structural pattern difference** only (convergence vs dispersion). We extend this into a **minimax information-theoretic envelope** (Math Thm 1): below $n\pi^2\Delta^2/(\sigma^2 d) \approx 9$, even the structural-pattern unlock cannot succeed. Path 1 implies "just detect the dispersion pattern" without a sample-complexity floor; our paper shows the floor binds empirically (9/9 retina holdouts below-floor confirm it). Our framework is therefore a strict superset.
- Path 1 implicitly treats protection as "once detected, you prevent it". Our pipeline separates detection (REAL) from protection (RareShield, Math Thm 5) formally. The separation matters because protection introduces its own identifiability questions (ELBO-preserving + convex-feasibility), which Path 1 doesn't touch.

## 6. What we still owe Criterion 3 (Scalability / Kang 4.9M)

This is the gap the human flagged as the immediate next priority. Our current evidence is:

- Synthetic pancreas factorial: Scanorama p=0.0099 on rare-truth motifs vs 0.32 on abundant.
- Pancreas real-data LossRate@10: Harmony 0.88 / Scanorama 0.98 / scVI 1.00 (above-floor regime).
- Retina real-data abundance series 9/9 at-floor per Thm 1 (below-floor regime).
- RareShield A/B on synth: AUPRC 0.556 → 0.867, ΔARI = 0.

**We have not yet run on Kang 2024 pan-cancer 4.9M.** Criterion 3 says the paper has limited practical value if we can't. The network blockers CompBio escalated earlier (CELLxGENE S3 checksum, figshare AWS WAF) are the binding constraint.

**Immediate team priority**: download Kang 2024 → 500k cancer-type subset first as a pilot → full 4.9M when feasible. Options that CompBio + ML + DS should try in parallel this cycle:

1. **cellxgene_census with payer context**: `cellxgene_census.open_soma(census_version="2024-07-01", context=SOMATileDBContext.create({"vfs.s3.request_payer": "requester"}))`.
2. **Manual figshare**: the dataset is listed on figshare with a known DOI; `wget` / `curl` with retries + resume, not the Python downloader. The figshare dataset ID for the Kang pan-cancer atlas is searchable.
3. **Zenodo mirror**: for the subset datasets; Baron/Muraro/Segerstolpe and HLCA have Zenodo mirrors; Kang may or may not.
4. **HuggingFace Datasets**: some pan-cancer atlases are re-hosted there.
5. **Human user manual stage**: place the H5AD at `/mnt/.../shared/data/kang_pancancer_immune.h5ad` (the path CompBio requested earlier).

Parallel to the download, ML + Math have already pushed the **hierarchical REAL** protocol — that is what we need to apply on Kang because the 4.9M atlas mixes cancer types and immune compartments that straddle the Theorem 1 floor. Full-atlas flat REAL will show the overcorrection pattern; hierarchical REAL (cancer-type stratified, then immune-compartment stratified) is expected to recover populations 100× rarer than the flat floor per Math Thm S1-hierarchical.

## 7. Concrete deliverables for the team this cycle

1. **CompBio**: Kang download via one of the five options above; pilot on one cancer-type slice (say LUAD or MEL) first; emit per-dataset H5AD with `adata.var['gene_symbol']` populated.
2. **DS**: Fig 4 real-data refresh as soon as a Kang slice lands; Fig 2 panels stay as-is.
3. **ML**: run RareShield on the Kang slice; report ΔARI on broad immune labels.
4. **Immunology + Tumor Biology**: 24h SLA on motif interpretation when CompBio ships Kang motif lists.
5. **Math**: empirical validation of Thm 1 + Thm S1-hierarchical on Kang — the paper's Criterion-3 scalability anchor.
6. **Philosophy**: audit any Kang claim against the 8-trap checklist + lint:epistemic before splice.
7. **Biological Science (me)**: splice Path 1's four convergent points into manuscript §01, §02, and §08 per §4 above; coordinate Kang arrival → v3.0 PDF.

## 8. v2.7 manuscript additions derived from Path 1

I'll implement the four additions from §4 in a separate v2.7 commit — small deltas, no new blockers:

- §01: add two sentences on the Triple-Absence self-closing loop.
- §02 Intuition: add one plain-language sentence "normal batch correction produces convergence; unannotated-rare-subpopulation erasure produces dispersion, and REAL's four components ($P_d, P_t, P_n, B$) measure exactly that distinction".
- §08 Limitations: add "circular-reasoning" language to the synth-data caveat.
- Abstract (optional): slight punchier phrasing per §4.4.

These are the content updates that do not depend on Kang arriving; they ship in v2.7 this cycle.
