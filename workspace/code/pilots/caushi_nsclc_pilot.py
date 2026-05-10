"""Caushi 2021 NSCLC × TPEX holdout pilot — CLONOTYPE-VALIDATED protocol.

STATUS: STUB. The Caushi 2021 NSCLC atlas (GSE176021, ~500k T cells across
~50 patients with paired TCR sequencing) has not yet been staged onto the
GPU server (4.3 GB tar files, manual download pending).

When data lands, this script implements the **third** (clonotype-anchored)
holdout modality per the Philosophy pre-audit checklist:

1. TCR clonotype labels are the validation channel (NOT MajorCluster). Each
   cell carries a TCR β/α chain string from the paired 10x TCR-seq output.
   "Same clonotype" = identical CDR3β sequence. TPEX expansion is typically
   measured by clonotype-level expansion from pre-treatment to on-treatment.

2. TRBV / TRAV / TRBJ / TRAJ gene-family members MUST be logged for HVG
   inclusion check. If any TRxV/J family gene lands in the top-2000 HVG,
   flag as Trap S0.1 (Philosophy register: "HVG-oracle-overlap"). If NONE
   land in HVG, the clonotype channel is genuinely orthogonal to the
   transcriptomic feature space — yielding the **strongest** non-circularity
   argument in the paper (stronger than Sade's P1/P2/P3 three-property guard,
   because TCR sequence identity is experimentally measured and lives outside
   the gene-expression matrix entirely).

3. Per-cell clonotype → motif assignment confusion matrix is reported, NOT
   per-cluster MajorCluster → motif confusion. This keeps the validation at
   the measurement-atomic level (individual TCR sequences per cell) rather
   than re-aggregating at cluster level.

Implementation sketch:

    import scirpy as ir

    def load_caushi_with_tcr(h5ad_path, tcr_clonotype_path):
        a = sc.read_h5ad(h5ad_path)
        # attach TCR clonotype obs column via ir.io.read_10x_vdj
        ir.io.read_10x_vdj(tcr_clonotype_path)  # adapted
        ir.pp.ir_dist(a)
        ir.tl.define_clonotypes(a, receptor_arms="all", dual_ir="primary_only")
        assert "clone_id" in a.obs
        return a

    def check_trxv_in_hvg(a, n_hvg=2000):
        sc.pp.highly_variable_genes(a, n_top_genes=n_hvg, batch_key="batch",
                                    flavor="seurat_v3")
        hvg = a.var.index[a.var["highly_variable"]]
        trxv = [g for g in hvg if g.startswith(("TRBV","TRAV","TRBJ","TRAJ",
                                                 "TRBC","TRAC","TRDV","TRDJ",
                                                 "TRGV","TRGJ"))]
        return trxv  # empty list = clean orthogonality

    def clonotype_motif_confusion(a, motif_key="leiden_preint_r3"):
        df = a.obs[["clone_id", motif_key]].copy()
        df = df[df["clone_id"].notna() & (df["clone_id"] != "no_IR")]
        return pd.crosstab(df["clone_id"], df[motif_key])

Output:
    shared/data/caushi_nsclc_tpex_holdout.h5ad
    shared/data/caushi_nsclc_tpex_holdout_TCR_sanity.json  # per the 3 checklist items

Pre-audit sign-off from Philosophy (agent-mozur8ub): YES if all three items
clean. Reviewer-tight positioning as "strongest non-circularity demonstration
in paper" iff HVG ∩ TCR-family genes = empty.
"""

from __future__ import annotations


def main() -> None:
    raise NotImplementedError(
        "Caushi 2021 NSCLC data not yet staged on /shared/data/. "
        "See Philosophy pre-audit checklist in this file's module docstring "
        "for the full protocol when the tar files land."
    )


if __name__ == "__main__":
    main()
