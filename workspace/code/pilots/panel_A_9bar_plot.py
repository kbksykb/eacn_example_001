"""Panel A appendix: 9-bar -log10(p_bh) method × scale visualization,
bars colored by ℓ_fine so the Scanorama-vs-Harmony signature (Math's ask) is
visible at a glance.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


def main() -> None:
    rows: list[dict] = []
    with open("workspace/results/locality_grid.csv") as f:
        for r in csv.DictReader(f):
            rows.append(r)

    dataset_order = ["cheng_paad", "cheng5", "cheng6"]
    method_order = ["harmony", "scanorama", "scvi"]
    dataset_label = {"cheng_paad": "Cheng PAAD\nn=2,853", "cheng5": "Cheng5\nn=20,341", "cheng6": "Cheng6\nn=49,271"}

    by_key = {(r["dataset"], r["method"]): r for r in rows}

    fig, ax = plt.subplots(figsize=(9, 5.2))

    bar_w = 0.24
    x_base = np.arange(len(dataset_order))

    cmap = plt.get_cmap("RdYlGn_r")  # low ℓ_fine (mixed) → red; high (preserved) → green-reversed? we want mixed=red
    # We want: high ℓ_fine = green (structure preserved), low = red (mixed). Use RdYlGn.
    cmap = plt.get_cmap("RdYlGn")

    method_offset = {"harmony": -bar_w, "scanorama": 0.0, "scvi": bar_w}
    method_hatch = {"harmony": "", "scanorama": "//", "scvi": "xx"}

    for method in method_order:
        xs = x_base + method_offset[method]
        heights = []
        colors = []
        for ds in dataset_order:
            r = by_key[(ds, method)]
            heights.append(float(r["neg_log10_pbh"]))
            colors.append(cmap(float(r["ell_fine_median"])))
        bars = ax.bar(xs, heights, width=bar_w, color=colors,
                      edgecolor="black", linewidth=0.8,
                      hatch=method_hatch[method], label=method)
        for b, h, ds in zip(bars, heights, dataset_order):
            r = by_key[(ds, method)]
            ax.text(b.get_x() + b.get_width() / 2, h + 0.06,
                    f"{h:.2f}", ha="center", va="bottom", fontsize=7)
            ax.text(b.get_x() + b.get_width() / 2, -0.18,
                    f"ℓf={float(r['ell_fine_median']):.2f}", ha="center", va="top", fontsize=6, color="#555")

    ax.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.8)
    ax.text(len(dataset_order) - 0.2, -np.log10(0.05) + 0.05, "α=0.05", fontsize=7, color="gray", ha="right")

    ax.set_xticks(x_base)
    ax.set_xticklabels([dataset_label[d] for d in dataset_order], fontsize=9)
    ax.set_ylabel("-log10(p_bh) on LAMP3+ mregDC rare-truth", fontsize=10)
    ax.set_title("Panel A appendix — 9-point method × scale grid\n"
                 "Bars colored by ℓ_fine (MajorCluster kNN locality; green=preserved, red=mixed)",
                 fontsize=10)
    ax.set_ylim(-0.45, 3.0)
    ax.grid(axis="y", alpha=0.3)

    legend_elts = [
        plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="black", hatch="", label="Harmony"),
        plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="black", hatch="//", label="Scanorama"),
        plt.Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="black", hatch="xx", label="scVI"),
    ]
    ax.legend(handles=legend_elts, loc="upper left", fontsize=9, title="Method")

    # Colorbar for ℓ_fine
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0.5, vmax=1.0))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("ℓ_fine (MajorCluster locality)", fontsize=9)

    fig.tight_layout()
    out = Path("workspace/results/panel_A_appendix_9bar.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
