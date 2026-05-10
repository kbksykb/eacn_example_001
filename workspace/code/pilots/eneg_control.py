"""
E_neg — batch-artefact negative control for REAL.

Philosophy's trap-7/trap-8 acceptance test: inject a cluster that LOOKS rare
(tight density mode) but is a pure batch-specific technical artefact, then
confirm REAL does NOT flag it as a loss.

Design:
- Base simulation as in synth_pilot (4 cell types, 3 batches, epsilon rare).
- Add E_neg cluster: in ONE batch only, shift 60 cells to a tight Gaussian
  cluster at a location that no other batch has cells near. This is a
  batch-specific technical artefact (doublets, ambient RNA bleed, batch-
  specific capture bias), NOT biology.
- Run integration (Harmony, Scanorama, scVI).
- Score REAL channels.
- The batch-artefact motif should:
    (a) NOT pass DS's cross-batch witness step (no other batch has matching signature).
    (b) Even if it slips through as a motif, the CoLM channel should NOT give a
        significant loss-probability because the pre- and post-integration
        densities at the artefact location differ only by smooth batch-effect
        deformation (which F admits).
"""

from __future__ import annotations

import argparse
import json
import logging
import pathlib
import time

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def simulate_with_eneg(seed: int, n_cells_per_batch: int = 2000, n_batches: int = 3,
                       eneg_n_cells: int = 60, eneg_shift_norm: float = 5.0):
    """Standard simulate + add a batch-0-only artefact cluster."""
    import anndata as ad
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    from workspace.code.datasets.synthetic import SimConfig

    # Base sim
    from workspace.code.datasets.synthetic import simulate
    a = simulate(SimConfig(seed=seed, n_cells_per_batch=n_cells_per_batch, n_batches=n_batches))

    # Add E_neg cluster in batch_0 only by selecting 60 cells from batch_0 and
    # injecting a tight per-cell shift in a housekeeping-gene direction.
    X = a.X.toarray() if hasattr(a.X, "toarray") else a.X
    is_b0 = (a.obs["batch"].astype(str) == "batch_0").to_numpy()
    b0_idx = np.where(is_b0)[0]
    pick = rng.choice(b0_idx, size=eneg_n_cells, replace=False)

    # Pick a specific housekeeping gene block (not a cell-type marker) and spike it high
    # — this creates a tight density mode in PCA space that doesn't match any biological type
    housekeeping_start = 400  # markers 0..399 are cell-type-specific, 400+ are shared/housekeeping
    spike_genes = np.arange(housekeeping_start, housekeeping_start + 50)
    X[np.ix_(pick, spike_genes)] += int(np.round(np.exp(eneg_shift_norm)))  # large integer Poisson bump

    a.X = sp.csr_matrix(X)
    # Mark the artefact membership
    is_eneg = np.zeros(a.n_obs, dtype=bool)
    is_eneg[pick] = True
    a.obs["is_eneg"] = is_eneg
    return a


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--methods", nargs="+", default=["harmony", "scanorama", "scvi"])
    p.add_argument("--out-name", default="synth_eneg")
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")

    from workspace.code.pilots.synth_pilot import METHODS, _pca_project, preprocess, over_cluster
    from workspace.code.real_channels import ChannelConfig, score_purity_and_procrustes
    from workspace.code.real_channels import ot_wrapper
    from workspace.code.real_channels.ot_channel import OTChannelConfig

    a = simulate_with_eneg(seed=args.seed)
    logger.info("simulated with E_neg: %d cells total, %d E_neg cells", a.n_obs, int(a.obs["is_eneg"].sum()))
    hidden = a.obs["__hidden_label__"].astype(str).to_numpy()

    a_pre = preprocess(a, batch_key="batch", seed=args.seed)
    candidates = over_cluster(a_pre, seed=args.seed, resolution=3.0)
    logger.info("motifs: %d", len(np.unique(candidates)))

    # Which motif corresponds to E_neg cluster?
    is_eneg = a.obs["is_eneg"].to_numpy()
    eneg_motif_ids = set()
    for cid in np.unique(candidates):
        m = candidates == cid
        if m.sum() > 0 and is_eneg[m].sum() / m.sum() > 0.5 and is_eneg[m].sum() >= 3:
            eneg_motif_ids.add(int(cid))
    logger.info("E_neg-dominated motifs: %s", sorted(eneg_motif_ids))

    # Rare-biological motif (epsilon)
    is_rare_bio = np.asarray([x == "epsilon" for x in hidden])
    bio_rare_ids = set()
    for cid in np.unique(candidates):
        m = candidates == cid
        if m.sum() >= 3 and is_rare_bio[m].sum() / max(m.sum(), 1) > 0.5 and is_rare_bio[m].sum() >= 3:
            bio_rare_ids.add(int(cid))
    logger.info("epsilon-dominated motifs: %s", sorted(bio_rare_ids))

    batch_int = pd.Categorical(a_pre.obs["batch"]).codes.astype(np.int32)
    cfg_ch = ChannelConfig(k=15, n_bootstrap=20, rare_abundance_threshold=0.10, seed=args.seed)
    ot_cfg = OTChannelConfig(n_permutations=100, device="cpu", k_neighbors=15, seed=args.seed)

    rows = []
    for method in args.methods:
        logger.info("=== %s ===", method)
        try:
            emb_int, meta = METHODS[method](a_pre, batch_key="batch", seed=args.seed)
        except Exception as exc:
            logger.exception("%s failed", method)
            continue
        d = emb_int.shape[1]
        pre_trunc = a_pre.obsm["X_uncorrected_pca"][:, :d]
        reports = score_purity_and_procrustes(
            emb_pre=pre_trunc.astype(np.float32),
            emb_post=emb_int.astype(np.float32),
            candidate_labels=candidates,
            batch_labels=batch_int,
            cfg=cfg_ch,
        )
        ot_res = ot_wrapper.score_two_sided(
            emb_pre=pre_trunc.astype(np.float32),
            emb_post=emb_int.astype(np.float32),
            candidate_labels=candidates,
            batch_labels=batch_int,
            cfg=ot_cfg,
        )
        ot_by_id = dict(zip(ot_res["candidate_ids"].tolist(), ot_res["p_values"].tolist()))

        for r in reports:
            rows.append({
                "method": method,
                "candidate_id": r.candidate_id,
                "abundance": r.abundance,
                "is_eneg_artefact": r.candidate_id in eneg_motif_ids,
                "is_biological_rare": r.candidate_id in bio_rare_ids,
                "loss_probability": r.loss_probability,
                "ot_p_value": ot_by_id.get(r.candidate_id, np.nan),
                "mknn_purity_post": r.purity_post,
            })

    df = pd.DataFrame(rows)
    print("\n=== E_neg negative-control result (should not flag artefact as loss) ===")
    print(df.to_string(index=False))

    out = pathlib.Path(__file__).resolve().parents[3] / "results" / f"{args.out_name}_neg_control.csv"
    df.to_csv(out, index=False)
    logger.info("wrote → %s", out)

    # Summary: for each method, is the E_neg artefact motif flagged?
    print("\n=== Summary ===")
    for m in args.methods:
        sub = df[df["method"] == m]
        eneg_sub = sub[sub["is_eneg_artefact"]]
        if len(eneg_sub) > 0:
            max_lp = float(eneg_sub["loss_probability"].max())
            min_ot = float(eneg_sub["ot_p_value"].min())
            print(f"{m:10s}  E_neg: max_loss_prob={max_lp:.3f}  min_ot_p={min_ot:.3f}  (should both be non-significant)")
        else:
            print(f"{m:10s}  E_neg: NO MOTIF CLAIMED THE ARTEFACT (correctly filtered by cross-batch witness)")


if __name__ == "__main__":
    main()
