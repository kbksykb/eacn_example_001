"""Plot the 6-point ℓ-locality regression (Math's falsifiable test)."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    rows: list[dict] = []
    with open("workspace/results/locality_grid.csv") as f:
        for r in csv.DictReader(f):
            rows.append(r)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    method_color = {"harmony": "#1f77b4", "scanorama": "#d62728", "scvi": "#2ca02c"}
    method_marker = {"harmony": "o", "scanorama": "s", "scvi": "^"}

    for ax, (ell_key, title) in zip(axes, [
        ("ell_coarse_times_nlog10p", "ℓ_coarse · -log10(p_bh)"),
        ("ell_fine_times_nlog10p", "ℓ_fine · -log10(p_bh)"),
    ]):
        for r in rows:
            x = float(r["n_pi2_delta2"])
            y = float(r[ell_key])
            m = r["method"]
            ax.scatter(x, y, color=method_color.get(m, "gray"),
                       marker=method_marker.get(m, "x"), s=110,
                       edgecolor="black", linewidth=0.7,
                       label=f"{m} × {r['dataset']}")
            ax.annotate(r["dataset"].replace("cheng_", "").replace("cheng", ""),
                        (x, y), xytext=(4, 4), textcoords="offset points", fontsize=7)
        ax.set_xscale("log")
        ax.set_xlabel("n·π²·Δ²  (log scale)")
        ax.set_ylabel(title)
        ax.set_title(title)
        ax.grid(alpha=0.3)

    # Single legend for the pair
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#1f77b4",
                   markeredgecolor="black", markersize=10, label="Harmony"),
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor="#d62728",
                   markeredgecolor="black", markersize=10, label="Scanorama"),
    ]
    axes[1].legend(handles=handles, loc="lower right", fontsize=9)

    fig.suptitle("Math's compartment-locality prediction (6-point grid). "
                 "Monotone collapse predicted if ℓ is the unifying variable.",
                 fontsize=10)
    fig.tight_layout()
    out = Path("workspace/results/locality_grid_6point.png")
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
