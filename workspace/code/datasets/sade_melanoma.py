"""Sade-Feldman 2018 melanoma anti-PD1 pre/post H5AD loader.

GSE120575. 16,291 single cells × 55,736 genes across 48 tumor samples from 32 patients.
Row 0 = cell_id; Row 1 = condition label (Pre_P1/Post_P1 etc.); Rows 2+ = gene_name + per-cell TPM.
Patient_ID file has the response label per patient.
"""

from __future__ import annotations

import gzip
import pathlib

import numpy as np
import pandas as pd


def main():
    src = pathlib.Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/data/sade2018")
    expr_path = src / "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"

    print("Parsing expression header rows...")
    with gzip.open(expr_path, "rt") as fh:
        cell_ids = fh.readline().strip().split("\t")  # row 0
        conditions = fh.readline().strip().split("\t")  # row 1 (Pre_P1, Post_P1, etc.)
        gene_names = []
        data_rows = []
        for line in fh:
            parts = line.strip().split("\t")
            gene_names.append(parts[0])
            data_rows.append([float(x) for x in parts[1:]])
    X = np.asarray(data_rows, dtype=np.float32).T  # cells × genes
    print(f"  cells × genes = {X.shape}")
    print(f"  first conditions: {conditions[:5]}")

    # Parse patient_ID file for response + sample metadata
    pid_path = src / "GSE120575_patient_ID_single_cells.txt.gz"
    pid_df = None
    try:
        with gzip.open(pid_path, "rt", encoding="utf-8", errors="replace") as fh:
            lines = fh.readlines()
        # Skip comment/empty/template lines; find the line starting with "title" (MIAMI template) or "Sample" etc.
        data_start = None
        for i, line in enumerate(lines):
            if line.startswith("title"):
                data_start = i
                break
        if data_start is not None:
            header = lines[data_start].strip().split("\t")
            rows = []
            for line in lines[data_start + 1 :]:
                parts = line.strip().split("\t")
                if not parts or not parts[0]:
                    continue
                rows.append(parts)
            max_cols = max(len(r) for r in rows) if rows else 0
            rows = [r + [""] * (max_cols - len(r)) for r in rows]
            if len(header) < max_cols:
                header = header + [f"col_{i}" for i in range(len(header), max_cols)]
            pid_df = pd.DataFrame(rows, columns=header[:max_cols])
            print(f"  patient_ID df: {pid_df.shape}, cols: {pid_df.columns.tolist()[:8]}")
    except Exception as exc:
        print(f"  patient_ID file read failed ({exc}); proceeding without response label")

    # Build obs
    obs = pd.DataFrame({"condition": conditions})
    obs.index = pd.Index(cell_ids, name="cell_id")
    # Parse condition into treatment + patient
    tp = obs["condition"].str.split("_", expand=True)
    obs["treatment"] = tp[0]  # Pre/Post
    obs["patient"] = tp[1]    # P1..P32
    obs["batch"] = obs["patient"].astype(str)
    obs["cell_type"] = "unknown"  # Sade-Feldman doesn't ship per-cell annotation in this file

    import anndata as ad
    import scipy.sparse as sp
    a = ad.AnnData(X=sp.csr_matrix(X))
    a.obs_names = pd.Index(cell_ids)
    a.var_names = pd.Index(gene_names)
    a.var["gene_symbol"] = a.var_names.astype(str).to_numpy()
    for col in obs.columns:
        a.obs[col] = obs[col].values

    out = src.parent / "sade_melanoma_cd8.h5ad"
    a.write_h5ad(out)
    print(f"\nWrote {out}")
    print(f"  shape {a.shape}")
    print(f"  treatment counts: {a.obs['treatment'].value_counts().to_dict()}")
    print(f"  patient count: {a.obs['patient'].nunique()}")


if __name__ == "__main__":
    main()
