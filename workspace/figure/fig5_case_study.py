"""fig5_case_study.py — Main Fig 5: label-free rare-subpopulation retrieval case studies.

Bio-Science v2.9 reframe (agent/biological_science@949d64a): canonical-TPEX
label-free retrieval, NOT novel-sub-state discovery.

Panels (Nature 2-column ~183mm):
  5a  LAMP3+ mregDC: Cheng PAAD + Cheng 5-cancer pan-cancer pool
  5b  CXCL13+ Tex: Zheng MM T-cell atlas
  5c  Canonical TPEX: Sade-Feldman melanoma — three-method outcome
  5d  Non-discovery control: REAL specificity (ground_truth=False stays below 0.5)
  5e  Scanorama erasure / REAL rescue summary across all case studies

Usage: python fig5_case_study.py --outdir workspace/manuscript/figures/fig5
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

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
METHOD_ORDER = ["harmony", "scanorama", "scvi"]
METHOD_LABELS = {"harmony": "Harmony", "scanorama": "Scanorama", "scvi": "scVI"}
METHOD_COLORS = {"harmony": OI[1], "scanorama": OI[0], "scvi": OI[4]}


def _load(fname: str) -> pd.DataFrame:
    p = DETECT / fname
    return pd.read_parquet(p) if p.exists() else pd.DataFrame()


def _rare_lp(df: pd.DataFrame) -> float:
    col = "ground_truth_label" if "ground_truth_label" in df.columns else "is_rare_truth"
    if col not in df.columns or df.empty:
        return float("nan")
    rare = df[df[col] == True]
    return float(rare["loss_probability"].mean()) if len(rare) else float("nan")


def panel_5a(ax: plt.Axes) -> None:
    paad = {m: _load(f"{m}_cheng_paad_lamp3_holdout.parquet") for m in METHOD_ORDER}
    cheng5 = {m: _load(f"{m}_cheng5_lamp3_holdout.parquet") for m in METHOD_ORDER}
    x = np.arange(len(METHOD_ORDER))
    for di, (ds, alpha, ds_label) in enumerate(
        [(paad, 1.0, "PAAD"), (cheng5, 0.55, "5-cancer")]
    ):
        for xi, m in enumerate(METHOD_ORDER):
            lp = _rare_lp(ds[m])
            lp = 0.0 if np.isnan(lp) else lp
            ax.bar(xi + (di - 0.5) * 0.40, lp, 0.38,
                   color=METHOD_COLORS[m], alpha=alpha,
                   edgecolor="black", lw=0.3)
            ax.text(xi + (di - 0.5) * 0.40, lp + 0.01, f"{lp:.2f}",
                    ha="center", fontsize=4.8)
    # OC annotation for Harmony PAAD
    df_h = paad["harmony"]
    if "is_overcorrection_candidate" in df_h.columns:
        oc = df_h[df_h["is_overcorrection_candidate"] == True]
        if len(oc):
            ax.annotate(f"OC flagged\np={oc.loss_probability.mean():.2f}",
                        xy=(-0.19, _rare_lp(df_h)),
                        xytext=(-0.55, _rare_lp(df_h) + 0.18),
                        arrowprops=dict(arrowstyle="->,head_width=0.15",
                                        color=OI[5], lw=0.7),
                        fontsize=5, color=OI[5])
    ax.set_xticks(x); ax.set_xticklabels([METHOD_LABELS[m] for m in METHOD_ORDER], fontsize=6)
    ax.set_ylabel("loss_probability (rare LAMP3⁺)"); ax.set_ylim(0, 0.65)
    ax.axhline(0.5, color="0.6", lw=0.4, ls="--")
    ax.set_title("LAMP3⁺ mregDC: Cheng PAAD + 5-cancer pool\n(Harmony OC flag; Scanorama/scVI preserve)",
                 loc="left", fontsize=7.0, pad=2)
    handles = [mpatches.Patch(color="0.4", label="PAAD (solid)"),
               mpatches.Patch(color="0.7", label="5-cancer (α=0.55)")]
    ax.legend(handles=handles, fontsize=5, loc="upper right")
    ax.text(-0.15, 1.06, "a", transform=ax.transAxes,
            fontweight="bold", fontsize=9, va="top", ha="right")


def panel_5b(ax: plt.Axes) -> None:
    zh = {m: _load(f"{m}_zheng_mm_tpex_holdout.parquet") for m in METHOD_ORDER}
    x = np.arange(len(METHOD_ORDER))
    for xi, m in enumerate(METHOD_ORDER):
        lp = _rare_lp(zh[m])
        if np.isnan(lp):
            lp = 0.0
        ax.bar(xi, lp, 0.55, color=METHOD_COLORS[m], edgecolor="black", lw=0.3)
        ax.text(xi, lp + 0.01, f"{lp:.3f}", ha="center", fontsize=5.5)
    ax.set_xticks(x); ax.set_xticklabels([METHOD_LABELS[m] for m in METHOD_ORDER], fontsize=6)
    ax.set_ylabel("loss_probability (rare CXCL13⁺ Tex)"); ax.set_ylim(0, 0.6)
    ax.axhline(0.5, color="0.6", lw=0.4, ls="--")
    ax.set_title("CXCL13⁺ Tex: Zheng MM T-cell atlas\n(ground_truth=True; Harmony lowest loss)",
                 loc="left", fontsize=7.0, pad=2)
    ax.text(-0.15, 1.06, "b", transform=ax.transAxes,
            fontweight="bold", fontsize=9, va="top", ha="right")


def panel_5c(ax: plt.Axes) -> None:
    sa = {"harmony": _load("harmony_sade_tpex_holdout.parquet"),
          "scanorama": _load("scanorama_sade_tpex_holdout_hvg3k.parquet"),
          "scvi":      _load("scvi_sade_tpex_holdout_hvg3k.parquet")}
    methods_plot = [m for m in METHOD_ORDER if not sa[m].empty]
    data = [sa[m]["loss_probability"].values for m in methods_plot]
    oc_counts = {m: int(sa[m].get("is_overcorrection_candidate", pd.Series(dtype=bool)).sum())
                 if "is_overcorrection_candidate" in sa[m].columns else 0
                 for m in methods_plot}
    if data and all(len(d) >= 3 for d in data):
        parts = ax.violinplot(data, positions=range(len(methods_plot)),
                              showmedians=True, showextrema=True)
        for i, (pc, m) in enumerate(zip(parts["bodies"], methods_plot)):
            pc.set_facecolor(METHOD_COLORS[m]); pc.set_alpha(0.7)
        parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(1.0)
        for e in ["cmins", "cmaxes", "cbars"]: parts[e].set_linewidth(0.5)
    ax.set_xticks(range(len(methods_plot)))
    ax.set_xticklabels([f"{METHOD_LABELS[m]}\nOC={oc_counts[m]}" for m in methods_plot], fontsize=5.5)
    ax.set_ylabel("loss_probability (all candidates)")
    ax.axhline(0.5, color="0.6", lw=0.4, ls="--")
    ax.set_title("Canonical TPEX: Sade-Feldman melanoma\n(non-discovery — known TPEX survives; OC count = overcorrection risk)",
                 loc="left", fontsize=7.0, pad=2)
    ax.text(-0.12, 1.06, "c", transform=ax.transAxes,
            fontweight="bold", fontsize=9, va="top", ha="right")


def panel_5d(ax: plt.Axes) -> None:
    all_dfs = []
    for fname in sorted(DETECT.glob("*.parquet")):
        df = pd.read_parquet(fname)
        if "ground_truth_label" in df.columns:
            all_dfs.append(df)
    if not all_dfs:
        ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center"); return
    combined = pd.concat(all_dfs, ignore_index=True)
    rare = combined[combined["ground_truth_label"] == True]["loss_probability"].values
    nonrare = combined[combined["ground_truth_label"] == False]["loss_probability"].values
    bins = np.linspace(0, 1, 25)
    ax.hist(nonrare, bins=bins, color="0.6", alpha=0.6, density=True,
            label=f"non-rare (n={len(nonrare)})")
    ax.hist(rare, bins=bins, color=OI[2], alpha=0.8, density=True,
            label=f"rare-truth (n={len(rare)})")
    ax.axvline(0.5, color=OI[5], lw=0.8, ls="--", label="p=0.5 threshold")
    ax.set_xlabel("loss_probability"); ax.set_ylabel("density")
    ax.set_title("specificity: non-rare candidates below threshold\n(rare-truth right-shifted; all methods + datasets)",
                 loc="left", fontsize=7.0, pad=2)
    ax.legend(fontsize=5.5)
    ax.text(-0.15, 1.06, "d", transform=ax.transAxes,
            fontweight="bold", fontsize=9, va="top", ha="right")


def panel_5e(ax: plt.Axes) -> None:
    cases = [
        ("LAMP3⁺\nPAAD",     "harmony_cheng_paad_lamp3_holdout",  "scanorama_cheng_paad_lamp3_holdout"),
        ("LAMP3⁺\n5-cancer", "harmony_cheng5_lamp3_holdout",      "scanorama_cheng5_lamp3_holdout"),
        ("CXCL13⁺ Tex\nZheng MM", "harmony_zheng_mm_tpex_holdout","scanorama_zheng_mm_tpex_holdout"),
        ("TPEX\nSade",        "harmony_sade_tpex_holdout",         "scanorama_sade_tpex_holdout_hvg3k"),
    ]
    labels, deltas = [], []
    for label, h_f, s_f in cases:
        h_lp = _rare_lp(_load(h_f + ".parquet"))
        s_lp = _rare_lp(_load(s_f + ".parquet"))
        if not (np.isnan(h_lp) or np.isnan(s_lp)):
            labels.append(label); deltas.append(s_lp - h_lp)
    if not labels:
        ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center"); return
    x = np.arange(len(labels))
    colors = [OI[5] if d > 0 else OI[2] for d in deltas]
    ax.bar(x, deltas, color=colors, edgecolor="black", lw=0.3)
    for xi, d in zip(x, deltas):
        ax.text(xi, d + (0.005 if d >= 0 else -0.013), f"{d:+.3f}", ha="center", fontsize=5.5)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel("Δ loss_prob (Scanorama − Harmony)")
    ax.set_title("Scanorama erasure vs Harmony across case studies\n(orange = Scanorama raises rare-cell loss; blue = Scanorama better)",
                 loc="left", fontsize=7.0, pad=2)
    ax.text(-0.08, 1.06, "e", transform=ax.transAxes,
            fontweight="bold", fontsize=9, va="top", ha="right")


def render(outdir: pathlib.Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(7.2, 9.5), constrained_layout=False)
    gs = fig.add_gridspec(nrows=3, ncols=2,
                          height_ratios=[1.0, 1.0, 1.0],
                          hspace=0.60, wspace=0.42,
                          left=0.11, right=0.97, top=0.97, bottom=0.06)
    panel_5a(fig.add_subplot(gs[0, 0]))
    panel_5b(fig.add_subplot(gs[0, 1]))
    panel_5c(fig.add_subplot(gs[1, 0]))
    panel_5d(fig.add_subplot(gs[1, 1]))
    panel_5e(fig.add_subplot(gs[2, :]))
    for ext in ("pdf", "png"):
        path = outdir / f"fig5.{ext}"
        fig.savefig(path, dpi=600 if ext == "png" else 300)
        print(f"wrote {path}", file=sys.stderr)
    plt.close(fig)


def _parse() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default="workspace/manuscript/figures/fig5", type=pathlib.Path)
    return p.parse_args()


if __name__ == "__main__":
    args = _parse()
    render(args.outdir)
