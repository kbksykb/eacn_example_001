"""ℓ-locality regression for Math's cross-method-compartment-amplification test.

For each (dataset, method) point, compute BOTH:
    ℓ_coarse := kNN locality at lineage resolution (DC / Mac / Mono / Mast / Neut)
    ℓ_fine   := kNN locality at MajorCluster resolution (M05_cDC3_LAMP3 etc.)

Math's prediction: ℓ × -log10(p_bh) should be monotone in n·π²·Δ² across ALL
(method, scale) points when the compartment-amplification mechanism is operative.
Hypothesis is that ℓ_fine is the discriminating metric — coarse DC/Mac/Mono
lineage separations are preserved by both methods, but fine intra-lineage
motif structure may dissolve differently under Harmony vs Scanorama.
"""

from __future__ import annotations

import csv
from pathlib import Path

import anndata as ad
import numpy as np
from sklearn.neighbors import NearestNeighbors

SHARED = Path("/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared")

GRID = [
    ("cheng_paad", "cheng_paad_lamp3_holdout.h5ad", "cheng_paad_myeloid.h5ad", 2853, 0.0107, 0.8, 1.0000, "harmony"),
    ("cheng5", "cheng5_lamp3_holdout.h5ad", "cheng_pancancer5_myeloid.h5ad", 20341, 0.0150, 0.8, 0.0130, "harmony"),
    ("cheng6", "cheng6_lamp3_holdout.h5ad", "cheng_pancancer_myeloid.h5ad", 49271, 0.0115, 0.8, 0.0041, "harmony"),
    ("cheng_paad", "cheng_paad_lamp3_holdout.h5ad", "cheng_paad_myeloid.h5ad", 2853, 0.0107, 0.8, 0.005, "scanorama"),
    ("cheng5", "cheng5_lamp3_holdout.h5ad", "cheng_pancancer5_myeloid.h5ad", 20341, 0.0150, 0.8, 0.008, "scanorama"),
    ("cheng6", "cheng6_lamp3_holdout.h5ad", "cheng_pancancer_myeloid.h5ad", 49271, 0.0115, 0.8, 0.0039, "scanorama"),
]


def coarse_lineage(label: str) -> str:
    if not isinstance(label, str):
        return "NA"
    if "cDC" in label or "LAMP3" in label:
        return "DC"
    if "Macro" in label:
        return "Mac"
    if "Mono" in label:
        return "Mono"
    if "Mast" in label:
        return "Mast"
    if "Neut" in label:
        return "Neut"
    return label.split("_")[0] if "_" in label else label


def ell_stats(X: np.ndarray, label: np.ndarray, k: int = 30) -> tuple[float, float, float]:
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm="auto", n_jobs=-1).fit(X)
    _, idx = nn.kneighbors(X)
    idx = idx[:, 1:]
    same = (label[idx] == label[:, None]).astype(np.float32)
    per_cell = same.mean(axis=1)
    return float(np.median(per_cell)), float(np.percentile(per_cell, 25)), float(np.percentile(per_cell, 75))


def main() -> None:
    out_rows: list[dict] = []
    lineage_cache: dict[str, tuple] = {}

    for dset, holdout, myel, n, pi, delta, p_bh, method in GRID:
        if dset not in lineage_cache:
            print(f"[load] lineage from {myel}")
            m = ad.read_h5ad(SHARED / "data" / myel, backed="r")
            major = m.obs["MajorCluster"].astype(str).to_numpy()
            coarse = np.array([coarse_lineage(x) for x in major])
            lineage_cache[dset] = (m.obs_names.to_numpy(), major, coarse)

        obs_m, major_m, coarse_m = lineage_cache[dset]

        int_path = SHARED / "integrations" / method / holdout
        print(f"[load] {int_path.name} ({method})")
        a = ad.read_h5ad(int_path)
        if "X_integrated" not in a.obsm:
            print(f"  MISSING X_integrated; skipping")
            continue

        obs_i = a.obs_names.to_numpy()
        map_idx = {b: i for i, b in enumerate(obs_m)}
        pos = np.array([map_idx.get(b, -1) for b in obs_i])
        mask = pos >= 0
        major_aligned = np.where(mask, major_m[pos.clip(0)], "NA")
        coarse_aligned = np.where(mask, coarse_m[pos.clip(0)], "NA")

        X = a.obsm["X_integrated"]
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = np.asarray(X, dtype=np.float32)

        ellc_med, ellc_q1, ellc_q3 = ell_stats(X, coarse_aligned, k=30)
        ellf_med, ellf_q1, ellf_q3 = ell_stats(X, major_aligned, k=30)

        nlog = float(-np.log10(max(p_bh, 1e-6)))
        n_pi2_d2 = float(n * pi * pi * delta * delta)

        print(f"  {dset} × {method}: ℓ_coarse={ellc_med:.3f}  ℓ_fine={ellf_med:.3f}  -log10(p)={nlog:.2f}  nπ²Δ²={n_pi2_d2:.2f}")

        out_rows.append(dict(
            dataset=dset, method=method, n=n, pi=pi, delta=delta,
            p_bh=p_bh, neg_log10_pbh=round(nlog, 4),
            ell_coarse_median=round(ellc_med, 4), ell_coarse_q1=round(ellc_q1, 4), ell_coarse_q3=round(ellc_q3, 4),
            ell_fine_median=round(ellf_med, 4), ell_fine_q1=round(ellf_q1, 4), ell_fine_q3=round(ellf_q3, 4),
            ell_coarse_times_nlog10p=round(ellc_med * nlog, 4),
            ell_fine_times_nlog10p=round(ellf_med * nlog, 4),
            n_pi2_delta2=round(n_pi2_d2, 4),
            rare_type="LAMP3_mregDC",
        ))

    outfp = Path("workspace/results/locality_grid.csv")
    outfp.parent.mkdir(parents=True, exist_ok=True)
    with outfp.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"wrote {outfp} with {len(out_rows)} rows")


if __name__ == "__main__":
    main()
