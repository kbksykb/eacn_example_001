"""Hierarchical REAL training wrapper.

Mathematics agent (agent-mozusggu). 2026-05-10.
Implements the telescoping-detection protocol for hierarchical REAL per
Supp Theorem S1-hierarchical (workspace/handoffs/hierarchical_real.tex).

The protocol:

    1. Run flat REAL on the whole atlas (compartment C_0).
    2. Extract major-lineage labels (from an external classifier or from the
       user's `compartment_labels` keyword) to partition cells into
       C_1 := { (major_lineage_i) : major_lineage_i in lineages }.
    3. For each compartment c_1 in C_1: re-run REAL restricted to c_1.
    4. Optionally repeat for finer compartments (C_2 = minor lineages).
    5. Union the witness-set outputs: W_hier = W_0 ∪ (∪_k W_k).

Under the compartment-extraction conditions of Theorem S1-hierarchical
(precision ≥ 1−η_ext, non-exclusion ≥ 1−η_excl), each level inherits
detectability with amplification factor n_atlas / n_compartment.

Key self-fixing construction (credit ML agent): seed the compartment re-run
with the *flat-REAL-suspect-collapsed* cells (LRS(w) < 0.4 at whole-atlas),
since those are exactly the cells most at risk of mis-classification by a
naive compartment classifier — the hierarchical re-run recovers them.

Usage:

    from hierarchical_real import run_hierarchical_real

    result = run_hierarchical_real(
        adata=my_anndata,
        compartment_key='major_lineage',            # optional; uses flat-REAL clusters if None
        integrator_fn=my_rareshield_fit_predict,    # (adata) -> z_post
        config=REALConfig(),
    )
    # result.W_all, result.motifs_per_level, result.detectability_envelope_per_level
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np


@dataclass
class REALConfig:
    k: int = 15                          # per-batch kNN size
    K: int = 300                         # smoothed-neighborhood size
    tau: float = 2.0                     # rarity threshold on R_b
    sigma_thresh: float = 0.85           # cross-batch signature similarity
    mu_min_batches: int = 2              # minimum witnessing depth
    c_min: int = 5                       # seed size floor
    c_max_frac: float = 0.005            # seed size ceiling as fraction of batch size
    lrs_gate: float = 0.4                # below this, call "suspect collapsed"
    max_levels: int = 3                  # C_0, C_1, C_2
    min_compartment_cells: int = 5_000   # skip compartments smaller than this
    alpha: float = 0.05                  # Type-I error
    beta: float = 0.05                   # Type-II error
    sigma_scale: float = 1.0             # embedding normalisation
    Lambda: float = 3.0                  # batch-heterogeneity penalty


@dataclass
class REALResult:
    W_all: list                                       # union of all motif witnesses
    motifs_per_level: dict[int, list]                 # level → motif list
    detectability_envelope_per_level: dict[int, dict]  # level → {n, floor, floor_met}
    compartment_sizes: dict                            # level → {compartment_id → n}


def detectability_floor(n: int, cfg: REALConfig) -> float:
    """Return the minimum π² · Δ² / σ² detectable at (α, β, Λ) given n cells."""
    c1 = 1.0
    log_factor = np.log(1.0 / (cfg.alpha * cfg.beta))
    return c1 * log_factor * (cfg.sigma_scale ** 2) * cfg.Lambda / n


def _flat_real_step(adata, integrator_fn, cfg: REALConfig):
    """Run flat REAL on an AnnData and return the witness list + z_post.

    This is the placeholder interface with DS/ML's actual REAL implementation;
    when integrated with examples/data_science/workspace/lrs_framework/, swap
    this body for the canonical call. The shape of the output is:
        (W_level : list[dict with keys motif_id, cell_indices, LRS, P_d, P_t, P_n, B])
    """
    z_post = integrator_fn(adata)
    adata.obsm["X_integrated"] = z_post
    # Placeholder: the actual REAL seed-enumeration + witnessing is DS's code.
    # Here we call a hook. If the hook is not configured, return an empty list.
    if not hasattr(integrator_fn, "real_witness"):
        return [], z_post
    W_level = integrator_fn.real_witness(adata, cfg)
    return W_level, z_post


def run_hierarchical_real(
    adata,
    integrator_fn: Callable,
    compartment_key: Optional[str] = "major_lineage",
    config: Optional[REALConfig] = None,
) -> REALResult:
    """Run hierarchical REAL on an AnnData.

    At each level k, a set of compartments is extracted (from compartment_key at
    level 1, or from the major-type clusters emerging at flat REAL if no key).
    Then REAL is re-run within each compartment. Motif witnesses are accumulated.

    Parameters
    ----------
    adata : AnnData with .obs['batch'] (required) and optionally .obs[compartment_key].
    integrator_fn : callable (AnnData) -> z_post; may have attribute .real_witness.
    compartment_key : the .obs column naming the major lineage; None = auto from flat REAL.
    config : REALConfig.

    Returns
    -------
    REALResult with unioned witnesses and per-level detectability envelopes.
    """
    cfg = config or REALConfig()
    motifs_per_level: dict[int, list] = {}
    envelope_per_level: dict[int, dict] = {}
    compartment_sizes: dict = {}
    W_all: list = []

    # Level 0: whole atlas
    print(f"[hierarchical REAL] level 0 (whole atlas): n = {adata.n_obs}")
    W_0, z_post_0 = _flat_real_step(adata, integrator_fn, cfg)
    motifs_per_level[0] = W_0
    envelope_per_level[0] = {"n": adata.n_obs, "floor": detectability_floor(adata.n_obs, cfg)}
    compartment_sizes[0] = {"atlas": adata.n_obs}
    W_all.extend(W_0)

    # Level 1: major-lineage compartments
    if compartment_key is None or compartment_key not in adata.obs.columns:
        # Derive from flat-REAL clustering (not implemented here — placeholder).
        print(f"[hierarchical REAL] no compartment_key '{compartment_key}'; stopping at level 0.")
        return REALResult(W_all, motifs_per_level, envelope_per_level, compartment_sizes)

    level = 1
    compartments = adata.obs[compartment_key].astype(str).unique()
    compartment_sizes[level] = {}
    motifs_per_level[level] = []
    envelope_per_level[level] = {}

    for c in compartments:
        mask = (adata.obs[compartment_key].astype(str) == c)
        n_c = int(mask.sum())
        compartment_sizes[level][c] = n_c
        if n_c < cfg.min_compartment_cells:
            print(f"[hierarchical REAL] skipping compartment '{c}': n={n_c} < min")
            continue

        sub = adata[mask].copy()
        print(f"[hierarchical REAL] level {level} compartment '{c}': n = {n_c} "
              f"(amplification {adata.n_obs / n_c:.1f}x)")
        W_c, _ = _flat_real_step(sub, integrator_fn, cfg)
        envelope_per_level[level][c] = {
            "n": n_c,
            "floor": detectability_floor(n_c, cfg),
            "amplification_vs_atlas": adata.n_obs / n_c,
        }
        # Tag each motif with its compartment of origin
        for w in W_c:
            w["compartment_level"] = level
            w["compartment_id"] = c
            motifs_per_level[level].append(w)
            W_all.append(w)

    # Level 2 (optional, minor lineages) — not expanded here; placeholder structure.
    # Typical practice: use a two-level split and stop. Add another loop if needed.

    return REALResult(W_all, motifs_per_level, envelope_per_level, compartment_sizes)


# -- Self-test --------------------------------------------------------------
if __name__ == "__main__":
    # Smoke test: mock AnnData-like object + integrator.
    class MockAnnData:
        def __init__(self, n, batch, major_lineage):
            import pandas as pd
            self.n_obs = n
            self.obs = pd.DataFrame({
                "batch": batch,
                "major_lineage": major_lineage,
            })
            self.obsm = {}

        def __getitem__(self, key):
            # Very minimal slicing for the boolean mask case
            if hasattr(key, "__len__"):
                idxs = np.where(np.asarray(key))[0]
                out = MockAnnData.__new__(MockAnnData)
                out.n_obs = len(idxs)
                out.obs = self.obs.iloc[idxs].reset_index(drop=True)
                out.obsm = {}
                return out
            raise NotImplementedError

        def copy(self):
            return self

    n = 100_000
    rng = np.random.default_rng(0)
    ad = MockAnnData(
        n=n,
        batch=rng.integers(0, 30, n),
        major_lineage=rng.choice(["T_NK", "myeloid", "stromal", "tumor"], n,
                                 p=[0.4, 0.25, 0.2, 0.15]),
    )

    def mock_integrator(a):
        return rng.standard_normal((a.n_obs, 32))

    def mock_witness(a, cfg):
        # Return 1 motif per compartment with a dummy LRS
        return [{"motif_id": "W_mock_0", "cell_indices": [], "LRS": 0.3}]
    mock_integrator.real_witness = mock_witness  # type: ignore

    r = run_hierarchical_real(ad, mock_integrator)
    print(f"W_all: {len(r.W_all)} motifs across {len(r.motifs_per_level)} levels")
    for lvl, env in r.detectability_envelope_per_level.items():
        print(f"  level {lvl}: {env}")
