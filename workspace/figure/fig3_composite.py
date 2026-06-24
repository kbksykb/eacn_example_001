"""fig3_composite.py — Fig 3 composite: real + synthetic ablation panels.

Combines real-data evidence with synthetic placeholders for panels
awaiting RareShield-trained H5ADs. Clearly labels each panel's data source.

Panels:
  3a  AUPRC per method (real: 1-loss_prob as rare-detection score across
      all 12 datasets with ground_truth labels). Harmony baseline.
  3b  Precision at loss_prob>0.5 threshold vs recall (real data).
  3c  RareShield ablation waterfall: Baseline 0.556 → F=identity 0.476
      → Full 0.917 (from Comp-Bio seed=0 synth — best real-data proxy
      available until full RareShield training lands).
  3d  RareShield A/B real retina: ΔARI + Δloss_prob on overcorrection
      candidates (from fig3_rareshield_ab.py data).
  3e  Per-dataset AUPRC heatmap (methods × datasets, real data).
  3f  Runtime scaling (real manifest.csv data).

Usage: python fig3_composite.py --outdir workspace/manuscript/figures/fig3
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from sklearn.metrics import average_precision_score, precision_recall_curve
    _HAS_SKLEARN = True
except ImportError:
    _HAS_SKLEARN = False

try:
    from lrs_framework import viz_style  # type: ignore
    viz_style.apply()
    OI = viz_style.OKABE_ITO
except Exception:  # noqa: BLE001
    OI = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
          "#0072B2", "#D55E00", "#CC79A7", "#000000"]
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 600, "savefig.bbox": "tight",
        "font.family": "sans-serif", "font.size": 7,
        "axes.titlesize": 8, "axes.labelsize": 7, "axes.linewidth": 0.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.labelsize": 6, "ytick.labelsize": 6,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2, "ytick.major.size": 2,
        "legend.fontsize": 6, "legend.frameon": False,
        "lines.linewidth": 0.8, "lines.markersize": 3, "pdf.fonttype": 42,
    })
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=OI)

DETECT = pathlib.Path("/mnt/eacn_example_001/shared/detections/real")
MANIFEST = pathlib.Path("/mnt/eacn_example_001/shared/integrations/manifest.csv")

METHOD_ORDER = ["harmony", "scanorama", "scvi", "bbknn"]
METHOD_LABELS = {"harmony": "Harmony", "scanorama": "Scanorama",
                 "scvi": "scVI", "bbknn": "BBKNN"}
METHOD_COLORS = {"harmony": OI[1], "scanorama": OI[0], "scvi": OI[4], "bbknn": OI[6]}

DATASET_SHORT = {
    "retina_bc89_holdout": "Retina\nBC8/9",
    "retina_bc4_holdout":  "Retina\nBC4",
    "retina_bc2_holdout":  "Retina\nBC2",
    "retina_bc1a_holdout": "Retina\nBC1a",
    "retina4k_bc89_holdout": "Retina\n4k",
    "cheng5_lamp3_holdout":  "Cheng5\nLAMP3",
    "cheng6_lamp3_holdout":  "Cheng6\nLAMP3",
    "cheng_paad_lamp3_holdout": "PAAD\nLAMP3",
    "pancreas_synth":       "Pancreas\n(synth)",
    "zheng_mm_tpex_holdout": "Zheng\nTPEX",
    "sade_tpex_holdout":    "Sade\nTPEX",
    "sade_tpex_holdout_hvg3k": "Sade\nHVG3k",
}


def _load_all() -> pd.DataFrame:
    parquets = sorted(DETECT.glob("*.parquet"))
    dfs = [pd.read_parquet(p) for p in parquets]
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def panel_3a(ax: plt.Axes, all_df: pd.DataFrame) -> None:
    """Overall AUPRC per method (1 − loss_prob as rare-detection score)."""
    if all_df.empty or not _HAS_SKLEARN:
        ax.text(0.5, 0.5, "data unavailable", transform=ax.transAxes, ha="center")
        return

    methods = [m for m in METHOD_ORDER if m in all_df["method"].values]
    auprc_vals, n_vals = [], []
    for m in methods:
        sub = all_df[(all_df["method"] == m) & all_df["ground_truth_label"].notna()].copy()
        y = (sub["ground_truth_label"] == True).astype(int).values
        s = (1 - sub["loss_probability"]).values
        if y.sum() == 0 or len(sub) < 3:
            auprc_vals.append(float("nan")); n_vals.append(0); continue
        auprc_vals.append(average_precision_score(y, s))
        n_vals.append(int(y.sum()))

    x = np.arange(len(methods))
    bars = ax.bar(x, auprc_vals,
                  color=[METHOD_COLORS[m] for m in methods],
                  edgecolor="black", lw=0.3)
    for xi, (v, n) in enumerate(zip(auprc_vals, n_vals)):
        if not np.isnan(v):
            ax.text(xi, v + 0.005, f"{v:.3f}\n(n={n})", ha="center", fontsize=5.3)

    # Random baseline
    prevalence = (all_df["ground_truth_label"] == True).mean()
    ax.axhline(prevalence, color="0.5", lw=0.5, ls="--",
               label=f"random baseline ({prevalence:.3f})")
    ax.set_xticks(x)
    ax.set_xticklabels([METHOD_LABELS[m] for m in methods], fontsize=6)
    ax.set_ylabel("AUPRC (rare-detection)")
    ax.set_ylim(0, 0.15)
    ax.set_title("REAL detection AUPRC per method\n"
                 "(score = 1 − loss_prob; rare-truth labels as ground truth)",
                 loc="left", fontsize=7.2, pad=2)
    ax.legend(fontsize=5.5, loc="upper right")
    ax.text(0.5, 0.04, "REAL data — all 12 datasets",
            transform=ax.transAxes, ha="center", fontsize=5.5, color=OI[2])
    ax.text(-0.12, 1.06, "a", transform=ax.transAxes,
            fontweight="bold", fontsize=9, va="top", ha="right")


def panel_3b(ax: plt.Axes, all_df: pd.DataFrame) -> None:
    """Precision vs recall at loss_prob > 0.5 threshold, per method."""
    if all_df.empty:
        ax.text(0.5, 0.5, "data unavailable", transform=ax.transAxes, ha="center")
        return

    methods = [m for m in METHOD_ORDER[:3] if m in all_df["method"].values]
    for m in methods:
        sub = all_df[(all_df["method"] == m) & all_df["ground_truth_label"].notna()]
        y = (sub["ground_truth_label"] == True).astype(int).values
        s = (1 - sub["loss_probability"]).values
        if y.sum() == 0 or not _HAS_SKLEARN:
            continue
        prec, rec, thr = precision_recall_curve(y, s)
        ax.step(rec, prec, where="post", color=METHOD_COLORS[m],
                lw=1.0, label=METHOD_LABELS[m], alpha=0.85)
        # Mark the loss_prob=0.5 operating point
        op_thr = 0.5  # score threshold = 1 − 0.5 = 0.5
        idx = np.searchsorted(thr[::-1], op_thr)
        idx = len(thr) - idx
        if 0 <= idx < len(prec):
            ax.scatter(rec[idx], prec[idx], s=25, color=METHOD_COLORS[m],
                       edgecolor="black", lw=0.4, zorder=5)

    # Random baseline
    prevalence = (all_df["ground_truth_label"] == True).mean()
    ax.axhline(prevalence, color="0.5", lw=0.5, ls="--",
               label=f"random ({prevalence:.3f})")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("recall"); ax.set_ylabel("precision")
    ax.set_title("precision–recall curves (rare-truth detection)\n"
                 "(dots = loss_prob > 0.5 operating point)",
                 loc="left", fontsize=7.2, pad=2)
    ax.legend(fontsize=5.5, loc="upper right")
    ax.text(0.5, 0.04, "REAL data — all 12 datasets",
            transform=ax.transAxes, ha="center", fontsize=5.5, color=OI[2])
    ax.text(-0.12, 1.06, "b", transform=ax.transAxes,
            fontweight="bold", fontsize=9, va="top", ha="right")


def panel_3c(ax: plt.Axes) -> None:
    """RareShield ablation waterfall (Comp-Bio seed=0 synthetic)."""
    labels = ["Baseline\n(scVI vanilla)", "F = identity\n(gate OFF)", "Full RareShield"]
    values = [0.556, 0.476, 0.917]
    colors = [OI[1], OI[5], OI[2]]

    x = np.arange(len(labels))
    ax.bar(x, values, color=colors, edgecolor="black", lw=0.3)
    for xi, v in zip(x, values):
        ax.annotate(f"{v:.3f}", (xi, v + 0.02), ha="center",
                    fontsize=6.5, fontweight="bold")

    # Delta arrows
    from matplotlib.patches import FancyArrowPatch
    d1 = values[1] - values[0]
    d2 = values[2] - values[1]
    ax.annotate(f"Δ={d1:+.3f}", xy=(0.5, min(values[0], values[1]) - 0.07),
                ha="center", fontsize=5.8, color=OI[5])
    ax.annotate(f"Δ={d2:+.3f}", xy=(1.5, min(values[1], values[2]) - 0.07),
                ha="center", fontsize=5.8, color=OI[2])
    ax.add_patch(FancyArrowPatch(
        (0.08, values[0]), (0.92, values[1]),
        arrowstyle="->,head_length=3,head_width=2",
        color=OI[5], lw=0.8, mutation_scale=5))
    ax.add_patch(FancyArrowPatch(
        (1.08, values[1]), (1.92, values[2]),
        arrowstyle="->,head_length=3,head_width=2",
        color=OI[2], lw=0.8, mutation_scale=5))

    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=5.8)
    ax.set_ylabel("AUPRC (seed=0 synth)"); ax.set_ylim(0, 1.05)
    ax.set_title("RareShield ablation waterfall\n"
                 "(gate OFF hurts −0.08; Full RareShield +0.44 vs gate-OFF)",
                 loc="left", fontsize=7.2, pad=2)
    ax.text(0.5, 0.04, "Comp-Bio seed=0 synthetic — RareShield full training pending",
            transform=ax.transAxes, ha="center", fontsize=5.3, color=OI[5], style="italic")
    ax.text(-0.12, 1.06, "c", transform=ax.transAxes,
            fontweight="bold", fontsize=9, va="top", ha="right")


def panel_3d(ax: plt.Axes) -> None:
    """RareShield A/B on real retina (Comp-Bio JSONs)."""
    # From fig3_rareshield_ab.py data
    runs = {
        "Harmony":   {"delta_ARI": 0.0037,  "mean_lp_delta_oc": 0.0,    "n_oc": 0},
        "Scanorama": {"delta_ARI": -0.0136, "mean_lp_delta_oc": -0.2357, "n_oc": 4},
    }
    gate_ari = 0.02
    methods = list(runs.keys())
    fig = ax.get_figure()

    # Sub-axes via inset — use a simple two-bar approach instead
    dA = [runs[m]["delta_ARI"] for m in methods]
    colors_a = [OI[2] if v >= 0 else OI[5] for v in dA]
    ax.bar(range(len(methods)), dA, color=colors_a, edgecolor="black", lw=0.3, width=0.4)
    ax.axhline(+gate_ari, color="0.5", lw=0.5, ls="--")
    ax.axhline(-gate_ari, color="0.5", lw=0.5, ls="--")
    for xi, v in enumerate(dA):
        ax.annotate(f"{v:+.4f}", (xi, v + (0.003 if v >= 0 else -0.006)),
                    ha="center", fontsize=5.5)

    # Overlay Δloss_prob as text annotation
    for xi, m in enumerate(methods):
        lp = runs[m]["mean_lp_delta_oc"]
        n = runs[m]["n_oc"]
        ax.text(xi, -0.055,
                f"Δlp_oc={lp:+.3f}\n(n={n})",
                ha="center", fontsize=5, color=OI[4])

    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(methods, fontsize=6)
    ax.set_ylabel("ΔARI (RareShield − vanilla)")
    ax.set_ylim(-0.07, 0.06)
    ax.set_title("RareShield A/B: real retina\n"
                 "(|ΔARI| < 0.02 gate passed; Scanorama lowers OC loss_prob by 0.24)",
                 loc="left", fontsize=7.2, pad=2)
    ax.text(0.5, 0.04, "REAL retina data (Comp-Bio A/B JSON)",
            transform=ax.transAxes, ha="center", fontsize=5.5, color=OI[2])
    ax.text(-0.12, 1.06, "d", transform=ax.transAxes,
            fontweight="bold", fontsize=9, va="top", ha="right")


def panel_3e(ax: plt.Axes, all_df: pd.DataFrame) -> None:
    """Per-dataset AUPRC heatmap (methods × datasets)."""
    if all_df.empty or not _HAS_SKLEARN:
        ax.text(0.5, 0.5, "data unavailable", transform=ax.transAxes, ha="center")
        return

    methods = [m for m in METHOD_ORDER[:3] if m in all_df["method"].values]
    datasets = sorted(all_df[all_df["ground_truth_label"].notna()]["dataset"].unique())
    # Only datasets that have at least one rare-truth candidate
    datasets = [d for d in datasets
                if (all_df[(all_df["dataset"] == d) &
                            (all_df["ground_truth_label"] == True)]).shape[0] > 0]

    mat = np.full((len(methods), len(datasets)), np.nan)
    for i, m in enumerate(methods):
        for j, d in enumerate(datasets):
            sub = all_df[(all_df["method"] == m) & (all_df["dataset"] == d) &
                         all_df["ground_truth_label"].notna()]
            y = (sub["ground_truth_label"] == True).astype(int).values
            s = (1 - sub["loss_probability"]).values
            if y.sum() == 0 or len(sub) < 3:
                continue
            try:
                mat[i, j] = average_precision_score(y, s)
            except Exception:
                pass

    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.5)
    ax.set_yticks(range(len(methods)))
    ax.set_yticklabels([METHOD_LABELS[m] for m in methods], fontsize=6)
    ax.set_xticks(range(len(datasets)))
    ax.set_xticklabels([DATASET_SHORT.get(d, d) for d in datasets],
                       rotation=35, ha="right", fontsize=5.3)
    ax.set_title("AUPRC heatmap: methods × datasets (real data)\n"
                 "(yellow = low, red = high; NaN = no rare-truth candidates)",
                 loc="left", fontsize=7.2, pad=2)
    cbar = plt.gcf().colorbar(im, ax=ax, shrink=0.75, pad=0.02)
    cbar.ax.tick_params(labelsize=5.5)
    cbar.set_label("AUPRC", fontsize=6)
    for i in range(len(methods)):
        for j in range(len(datasets)):
            if not np.isnan(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                        fontsize=4.5, color="black")
    ax.text(-0.15, 1.06, "e", transform=ax.transAxes,
            fontweight="bold", fontsize=9, va="top", ha="right")


def panel_3f(ax: plt.Axes) -> None:
    """Runtime scaling from real manifest.csv."""
    if not MANIFEST.exists():
        ax.text(0.5, 0.5, "manifest not found", transform=ax.transAxes, ha="center")
        return
    mf = pd.read_csv(MANIFEST)
    methods_plot = ["harmony", "scanorama", "scvi"]
    for m in methods_plot:
        sub = mf[mf["method"] == m].sort_values("n_cells_used")
        if len(sub) < 2:
            continue
        ax.scatter(sub["n_cells_used"], sub["wall_clock_seconds"],
                   color=METHOD_COLORS[m], s=12, zorder=3,
                   label=METHOD_LABELS[m])
        x = sub["n_cells_used"].values
        y = sub["wall_clock_seconds"].values
        if len(x) >= 2 and np.all(y > 0):
            slope = np.polyfit(np.log10(x), np.log10(y), 1)
            x_fit = np.logspace(np.log10(x.min()), np.log10(x.max()), 50)
            ax.loglog(x_fit, 10 ** np.polyval(slope, np.log10(x_fit)),
                      color=METHOD_COLORS[m], lw=0.8, alpha=0.6)
    ax.set_xlabel("cells in dataset")
    ax.set_ylabel("wall-clock (s)")
    ax.set_title("runtime scaling (real measured values from manifest.csv)",
                 loc="left", fontsize=7.2, pad=2)
    ax.legend(fontsize=5.5, loc="upper left")
    ax.grid(True, which="both", ls=":", lw=0.3, alpha=0.5)
    ax.text(-0.16, 1.06, "f", transform=ax.transAxes,
            fontweight="bold", fontsize=9, va="top", ha="right")


def render(outdir: pathlib.Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    all_df = _load_all()
    # Normalise ground_truth_label
    if not all_df.empty and "ground_truth_label" not in all_df.columns and "is_rare_truth" in all_df.columns:
        all_df = all_df.rename(columns={"is_rare_truth": "ground_truth_label"})

    n_real = int((all_df["ground_truth_label"] == True).sum()) if not all_df.empty else 0
    print(f"Loaded {len(all_df)} rows, {n_real} rare-truth candidates", file=sys.stderr)

    fig = plt.figure(figsize=(7.2, 9.5), constrained_layout=False)
    gs = fig.add_gridspec(
        nrows=3, ncols=2,
        height_ratios=[1.0, 1.0, 1.1],
        hspace=0.60, wspace=0.40,
        left=0.10, right=0.97, top=0.97, bottom=0.07,
    )
    panel_3a(fig.add_subplot(gs[0, 0]), all_df)
    panel_3b(fig.add_subplot(gs[0, 1]), all_df)
    panel_3c(fig.add_subplot(gs[1, 0]))
    panel_3d(fig.add_subplot(gs[1, 1]))
    panel_3e(fig.add_subplot(gs[2, 0]), all_df)
    panel_3f(fig.add_subplot(gs[2, 1]))

    for ext in ("pdf", "png"):
        path = outdir / f"fig3_composite.{ext}"
        fig.savefig(path, dpi=600 if ext == "png" else 300)
        print(f"wrote {path}", file=sys.stderr)
    plt.close(fig)


def _parse() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default="workspace/manuscript/figures/fig3",
                   type=pathlib.Path)
    return p.parse_args()


if __name__ == "__main__":
    args = _parse()
    render(args.outdir)
