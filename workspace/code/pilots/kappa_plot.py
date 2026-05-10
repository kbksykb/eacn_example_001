"""Plot κ grid vs -log10(p_bh), bars per method × scale."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    rows = []
    with open("workspace/results/kappa_grid.csv") as f:
        for r in csv.DictReader(f):
            rows.append(r)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))

    method_color = {"harmony": "#1f77b4", "scanorama": "#d62728", "scvi": "#2ca02c"}
    method_marker = {"harmony": "o", "scanorama": "s", "scvi": "^"}
    dset_x = {"cheng_paad": 2853, "cheng5": 20341, "cheng6": 49271, "zheng_mm": 12274}

    # Left: κ vs n (log x)
    for m in ["harmony", "scanorama", "scvi"]:
        xs, ys, errs_low, errs_high = [], [], [], []
        for r in rows:
            if r["method"] != m:
                continue
            n = dset_x[r["dataset"]]
            xs.append(n)
            ys.append(float(r["kappa_median"]))
            errs_low.append(float(r["kappa_median"]) - float(r["kappa_q1"]))
            errs_high.append(float(r["kappa_q3"]) - float(r["kappa_median"]))
        order = np.argsort(xs)
        xs = np.array(xs)[order]
        ys = np.array(ys)[order]
        errs = np.vstack([np.array(errs_low)[order], np.array(errs_high)[order]])
        ax1.errorbar(xs, ys, yerr=errs, marker=method_marker[m], ms=9, lw=2,
                     color=method_color[m], label=m, capsize=4, alpha=0.85)
    ax1.set_xscale("log")
    ax1.set_xlabel("n_cells (log)")
    ax1.set_ylabel("κ (motif centroid pre/post displacement, median ± IQR)")
    ax1.set_title("κ is nearly flat across scale for each method")
    ax1.grid(alpha=0.3)
    ax1.legend()

    # Right: κ vs -log10(p_bh)
    for r in rows:
        m = r["method"]
        ax2.scatter(float(r["kappa_median"]), float(r["neg_log10_pbh"]),
                    color=method_color[m], marker=method_marker[m], s=110,
                    edgecolor="black", linewidth=0.7)
        ax2.annotate(r["dataset"].replace("cheng_", "").replace("cheng", ""),
                     (float(r["kappa_median"]), float(r["neg_log10_pbh"])),
                     xytext=(5, 5), textcoords="offset points", fontsize=7)
    ax2.set_xlabel("κ_median (per-method cross-batch redistribution)")
    ax2.set_ylabel("-log10(p_bh) on LAMP3+ rare-truth")
    ax2.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
    ax2.set_title("κ stratifies detection regimes: low κ=scale-mediated, high κ=method-intrinsic")
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    out = Path("workspace/results/kappa_grid_9point.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
