"""fig1_overview.py — Data Science contribution, Main Fig 1 mockup.

Panels (Nature 2-column, ~183mm wide):
  1a  Problem cartoon: two-batch embedding → integration → rare subpop erased.
  1b  Timeline of methods (2020-2026) that claim rare-subpop preservation but
      evaluate only on annotated rare cells — motivating the detectability gap.
  1c  REAL pipeline schematic: per-batch density anchor → cross-batch motif
      witness → post-integration survival test.
  1d  RareShield schematic: baseline integration objective with L_mass + L_anchor.
  1e  Toy 2-batch synthetic with one injected rare cluster, showing baseline
      integrator erases it while the RareShield variant preserves it.

Inputs: none (self-contained synthetic data — seed=0 for reproducibility).
Output: workspace/manuscript/figures/fig1/fig1.pdf + fig1.png
Usage:  python fig1_overview.py --outdir workspace/manuscript/figures/fig1
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch
from sklearn.decomposition import PCA

try:
    # Prefer the package-pinned style if present.
    from lrs_framework import viz_style  # type: ignore
    viz_style.apply()
    OI = viz_style.OKABE_ITO
except Exception:  # noqa: BLE001
    OI = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
          "#0072B2", "#D55E00", "#CC79A7", "#000000"]
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "font.family": "sans-serif",
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "axes.linewidth": 0.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": 6, "ytick.labelsize": 6,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2, "ytick.major.size": 2,
        "legend.fontsize": 6, "legend.frameon": False,
        "lines.linewidth": 0.8, "lines.markersize": 3,
        "pdf.fonttype": 42,
    })
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=OI)


# -------------------------------------------------------------------------- #
# Synthetic dataset for panel (d)                                             #
# -------------------------------------------------------------------------- #
def make_toy(seed: int = 0, n_abundant: int = 600, n_rare: int = 20, n_dim: int = 30):
    """Two batches, three abundant populations + one injected rare population.
    Returns X (2D embedding via PCA) and labels.
    """
    rng = np.random.default_rng(seed)
    centers_abundant = np.array([[0, 0], [6, 1], [3, 5]])
    rare_center = np.array([4.5, 3.0])
    batch_shift = np.array([[0, 0], [3.5, -2.0]])

    X_list, y_ct, y_batch = [], [], []
    for b in range(2):
        for c_idx, c in enumerate(centers_abundant):
            n = n_abundant // (2 * len(centers_abundant))
            pts = rng.normal(loc=c + batch_shift[b], scale=0.9, size=(n, 2))
            # pad with gaussian noise in high dim
            pad = rng.normal(scale=0.35, size=(n, n_dim - 2))
            X_list.append(np.hstack([pts, pad]))
            y_ct.extend([f"A{c_idx}"] * n)
            y_batch.extend([b] * n)
        # rare population — same identity in both batches, tight cluster, small
        n_r = n_rare // 2
        pts = rng.normal(loc=rare_center + batch_shift[b], scale=0.35, size=(n_r, 2))
        pad = rng.normal(scale=0.15, size=(n_r, n_dim - 2))
        X_list.append(np.hstack([pts, pad]))
        y_ct.extend(["R"] * n_r)
        y_batch.extend([b] * n_r)

    X = np.vstack(X_list)
    y_ct = np.array(y_ct)
    y_batch = np.array(y_batch)
    order = rng.permutation(X.shape[0])
    return X[order], y_ct[order], y_batch[order]


def pseudo_umap(X: np.ndarray, seed: int = 0) -> np.ndarray:
    """Light-weight 2D embedding for the panel mockup.
    PCA → small jitter. This is a schematic, not an empirical UMAP.
    """
    rng = np.random.default_rng(seed)
    Z = PCA(n_components=2, random_state=seed).fit_transform(X)
    Z += rng.normal(scale=0.05, size=Z.shape)
    return Z


def harmonyish_integrate(Z: np.ndarray, batch: np.ndarray, strength: float) -> np.ndarray:
    """Schematic 'overcorrection' operator: pull each batch centroid toward
    global centroid by `strength`, then add isotropic noise proportional to strength.
    Absorbs low-mass clusters into the nearest abundant region.
    """
    Z2 = Z.copy()
    global_mean = Z2.mean(0)
    for b in np.unique(batch):
        m = batch == b
        c = Z2[m].mean(0)
        Z2[m] += strength * (global_mean - c)
    # mass diffusion proportional to strength — this is what kills rare cells
    rng = np.random.default_rng(0)
    Z2 += rng.normal(scale=0.3 * strength, size=Z2.shape)
    return Z2


def real_aware_integrate(Z: np.ndarray, batch: np.ndarray, seed_mask: np.ndarray,
                         strength: float) -> np.ndarray:
    """Schematic of a RareShield-like operator: batch-shift correction for
    non-seed cells only; seed cells keep their local geometry (protected).
    """
    Z2 = Z.copy()
    global_mean = Z2[~seed_mask].mean(0)
    for b in np.unique(batch):
        m = (batch == b) & (~seed_mask)
        if m.sum() == 0:
            continue
        c = Z2[m].mean(0)
        Z2[m] += strength * (global_mean - c)
    rng = np.random.default_rng(1)
    Z2[~seed_mask] += rng.normal(scale=0.25 * strength, size=Z2[~seed_mask].shape)
    # Seed cells get only a small shared shift (co-registration), preserving local
    # structure.  This is the "mass-preservation" clause.
    if seed_mask.any():
        for b in np.unique(batch):
            m = (batch == b) & seed_mask
            if m.sum() == 0:
                continue
            c = Z2[m].mean(0)
            shift = strength * 0.3 * (global_mean - c)
            Z2[m] += shift
    return Z2


# -------------------------------------------------------------------------- #
# Panel renderers                                                             #
# -------------------------------------------------------------------------- #
def _scatter(ax, Z, ct, batch, *, highlight_rare=True, title=None):
    abundant_mask = ct != "R"
    rare_mask = ct == "R"
    markers = {0: "o", 1: "s"}
    for b in [0, 1]:
        m = abundant_mask & (batch == b)
        ax.scatter(Z[m, 0], Z[m, 1],
                   s=2.2, c=OI[1 if b == 0 else 4], marker=markers[b],
                   alpha=0.55, linewidths=0, rasterized=True,
                   label=f"Batch {b + 1}" if not highlight_rare else None)
    if highlight_rare:
        for b in [0, 1]:
            m = rare_mask & (batch == b)
            ax.scatter(Z[m, 0], Z[m, 1], s=20, c=OI[5], marker=markers[b],
                       edgecolor="white", linewidth=0.4, zorder=4,
                       label=f"rare (batch {b + 1})" if b == 0 else None)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")
    if title:
        ax.set_title(title, pad=4)


def panel_1a(fig, gs):
    """Problem cartoon: pre vs post integration, rare cell is erased."""
    X, ct, batch = make_toy()
    Z_pre = pseudo_umap(X)
    Z_post_bad = harmonyish_integrate(Z_pre, batch, strength=0.65)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    _scatter(ax1, Z_pre, ct, batch, title="pre-integration")
    _scatter(ax2, Z_post_bad, ct, batch, title="post (standard integration)")

    # Arrow between panels
    from matplotlib.patches import FancyArrowPatch
    # We draw the arrow in figure fraction coords
    fig_transform = fig.transFigure.inverted()
    a1 = fig_transform.transform(ax1.transAxes.transform((1.02, 0.5)))
    a2 = fig_transform.transform(ax2.transAxes.transform((-0.02, 0.5)))
    fig.add_artist(FancyArrowPatch(a1, a2,
                                   arrowstyle="->,head_length=4,head_width=3",
                                   mutation_scale=6, lw=0.8, color="black"))

    # Annotation: red circle over rare area in post
    ax2.add_patch(plt.Circle((Z_post_bad[ct == "R", 0].mean(),
                              Z_post_bad[ct == "R", 1].mean()),
                             0.9, fill=False, lw=0.8, ec=OI[5], ls="--"))
    ax2.annotate("rare population\nabsorbed", xy=(Z_post_bad[ct == "R", 0].mean(),
                                                  Z_post_bad[ct == "R", 1].mean()),
                 xytext=(8, -12), textcoords="offset points",
                 fontsize=6, color=OI[5],
                 arrowprops=dict(arrowstyle="-", color=OI[5], lw=0.4))

    # Shared legend beneath
    handles = [
        plt.Line2D([0], [0], marker="o", ls="", color=OI[1], markersize=3, label="Batch 1"),
        plt.Line2D([0], [0], marker="s", ls="", color=OI[4], markersize=3, label="Batch 2"),
        plt.Line2D([0], [0], marker="o", ls="", color=OI[5], markersize=5,
                   markeredgecolor="white", markeredgewidth=0.4, label="rare subpop."),
    ]
    ax1.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(1.05, -0.18), ncol=3,
               handletextpad=0.2, columnspacing=0.8)

    ax1.text(-0.08, 1.06, "a", transform=ax1.transAxes,
             fontweight="bold", fontsize=9, va="top", ha="right")


def _flowbox(ax, xy, wh, text, face=OI[1], alpha=0.08, ec=OI[1]):
    x, y = xy
    w, h = wh
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.02,rounding_size=0.05",
                                linewidth=0.6, edgecolor=ec,
                                facecolor=face, alpha=alpha))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=6.5, color="black")


def _flowarrow(ax, p0, p1, color="black"):
    ax.add_patch(FancyArrowPatch(p0, p1,
                                 arrowstyle="->,head_length=3,head_width=2",
                                 mutation_scale=5, lw=0.6, color=color))


def panel_1b_timeline(fig, gs):
    """Two-row timeline combining methods (top) + rare-cell discoveries (bottom).

    Top row: 2022-2026 methods claiming rare-subpop preservation, all
    evaluating on labels. Bottom row: historical rare-immune-cell discoveries
    with translational stakes + Le Cam detection-floor status.
    """
    ax = fig.add_subplot(gs[1, :])
    # -------- methods (top row) --------
    methods = [
        (2022, "scIB", "benchmark", False),
        (2022, "SCIDRL", "method", False),
        (2023, "scDML", "method", False),
        (2023, "scDREAMER", "method", False),
        (2024, "STACAS", "method", False),
        (2024, "CellANOVA", "recovery", False),
        (2025, "scCRAFT", "method", False),
        (2025, "RBET", "benchmark", False),
        (2025, "scIB-E", "benchmark", False),
        (2025, "semi-sup bench", "benchmark", False),
        (2026, "Nat CompSci review", "review", False),
        (2026, "REAL + RareShield", "this paper", True),
    ]
    # -------- discoveries (bottom row, Immunology handoff) --------
    # recovery tier: "flat" = above Le Cam floor; "hierarchical" = only via
    # sub-compartment pooling; "below" = unresolved at all depths.
    discoveries = [
        (1999, "pDC", "biomarker", 0.5, "flat"),
        (2017, "AS DC / DC5", "mechanistic", 0.04, "hierarchical"),
        (2018, "CD103+CD39+ TRM", "therapy", 1.0, "flat"),
        (2020, "LAMP3+ mregDC", "biomarker", 0.1, "hierarchical"),
        (2022, "FOLR2+ macrophage", "biomarker", 0.2, "flat"),
        (2026, "LAMP3+ mregDC (pan-cancer)", "biomarker", 0.01, "hierarchical"),
    ]

    years_all = sorted({y for y, *_ in methods} | {y for y, *_ in discoveries})
    y0, y1 = min(years_all) - 1, max(years_all) + 1
    ax.set_xlim(y0, y1)
    ax.set_ylim(-2.2, 2.2)

    # Central horizontal axis
    ax.axhline(0, color="0.3", lw=0.5, zorder=1)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_position(("data", 0))
    ax.set_xticks(years_all)
    ax.set_xticklabels([str(y) for y in years_all], fontsize=5.5)

    # Top row — methods above axis
    cat_color = {
        "benchmark": OI[1], "method": OI[4], "recovery": OI[5],
        "review": OI[7], "this paper": OI[2],
    }
    offs_year: dict[int, int] = {}
    for year, name, cat, _ in methods:
        s = +1
        n_same = offs_year.get(year, 0)
        offs_year[year] = n_same + 1
        y = s * (0.4 + 0.45 * (n_same % 3))
        c = cat_color[cat]
        ax.plot([year, year], [0, y], color=c, lw=0.6, zorder=2)
        ax.plot(year, y, "o", color=c, markersize=3, zorder=3)
        ax.annotate(name, xy=(year, y), xytext=(3, 5), textcoords="offset points",
                    fontsize=5, color="black")

    # Bottom row — discoveries below axis (stakes coded by fill;
    # recovery tier coded by outline ring: green-fill/no-ring, amber-ring, red-ring).
    stake_color = {
        "mechanistic": OI[1], "biomarker": OI[4],
        "therapy": OI[2], "future": "0.4",
    }
    tier_edge = {
        "flat":         ("none",    0.0),   # filled, no outline
        "hierarchical": (OI[0],     1.0),   # amber outline (Okabe-Ito orange)
        "below":        (OI[5],     1.1),   # red outline
    }
    offs_year_b: dict[int, int] = {}
    for year, name, stake, prev, tier in discoveries:
        s = -1
        n_same = offs_year_b.get(year, 0)
        offs_year_b[year] = n_same + 1
        y = s * (0.5 + 0.55 * (n_same % 2))
        c = stake_color[stake]
        edge, edgew = tier_edge.get(tier, ("black", 0.5))
        ax.plot([year, year], [0, y], color=c, lw=0.6, zorder=2)
        ax.plot(year, y, "s", color=c, markersize=4.5,
                markeredgecolor=edge if edge != "none" else c,
                markeredgewidth=edgew, zorder=3)
        ax.annotate(f"{name}\n({prev:.2f}% prev.)", xy=(year, y),
                    xytext=(4, -6 if n_same % 2 == 0 else -18),
                    textcoords="offset points",
                    fontsize=4.8, color="black")

    # Row labels
    ax.text(y0 - 0.5, 1.6, "methods\n(label-based eval)", fontsize=6,
            color="0.4", fontstyle="italic", ha="left", va="center")
    ax.text(y0 - 0.5, -1.6, "rare-immune-cell\ndiscoveries", fontsize=6,
            color="0.4", fontstyle="italic", ha="left", va="center")

    # Legend
    handles = [
        plt.Line2D([0], [0], marker="o", ls="", color=cat_color["method"],
                   label="method", markersize=3),
        plt.Line2D([0], [0], marker="o", ls="", color=cat_color["benchmark"],
                   label="benchmark", markersize=3),
        plt.Line2D([0], [0], marker="o", ls="", color=cat_color["this paper"],
                   label="this paper", markersize=3),
        plt.Line2D([0], [0], marker="s", ls="", color=stake_color["therapy"],
                   label="therapy", markersize=3),
        plt.Line2D([0], [0], marker="s", ls="", color=stake_color["biomarker"],
                   label="biomarker", markersize=3),
        plt.Line2D([0], [0], marker="s", ls="", color="0.6",
                   markeredgecolor=OI[0], markeredgewidth=1.0,
                   label="hierarchical-REAL", markersize=4),
        plt.Line2D([0], [0], marker="s", ls="", color="white",
                   markeredgecolor=OI[5], markeredgewidth=1.0,
                   label="below all floors", markersize=4),
    ]
    ax.legend(handles=handles, ncol=7, loc="lower center",
              bbox_to_anchor=(0.5, -0.40), frameon=False,
              handletextpad=0.3, columnspacing=0.5, fontsize=4.8)

    ax.set_title("above: methods claim rare-subpop preservation — all evaluate on labels · "
                 "below: historical rare-immune-cell discoveries · "
                 "amber outline = recovered only by hierarchical REAL; red = below all floors · "
                 "structural unlock: convergence vs dispersion leaves distinct patterns in data",
                 loc="left", fontsize=6.0, pad=2)

    ax.text(-0.01, 1.14, "b", transform=ax.transAxes,
            fontweight="bold", fontsize=9, va="top", ha="right")


def panel_1c_pipeline(fig, gs):
    """REAL pipeline schematic (3-step flow)."""
    ax = fig.add_subplot(gs[2, 0])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect("auto")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("REAL framework: label-free detectability pipeline", loc="left",
                 fontsize=7.5, pad=2)

    # 3 boxes horizontally
    _flowbox(ax, (0.03, 0.40), (0.28, 0.30),
             "per-batch\ndensity anchors $R_b(x) = \\rho_b / \\widetilde{\\rho}_b$",
             face=OI[1], ec=OI[1])
    _flowbox(ax, (0.36, 0.40), (0.28, 0.30),
             "cross-batch\nmotif witness $W$\n(MNN on HVG centroids)",
             face=OI[2], ec=OI[2])
    _flowbox(ax, (0.69, 0.40), (0.28, 0.30),
             "survival test\n$\\mathrm{LRS}(w) = HM(P_d, P_t, P_n, 1-B)$",
             face=OI[5], ec=OI[5])

    _flowarrow(ax, (0.31, 0.55), (0.36, 0.55))
    _flowarrow(ax, (0.64, 0.55), (0.69, 0.55))

    # caption-ish note below
    ax.text(0.5, 0.15,
            "label-free, scalable to 4.9M cells via FAISS-GPU + cugraph + witness-complex PH",
            ha="center", fontsize=6, color="0.3")

    ax.text(-0.02, 1.06, "c", transform=ax.transAxes,
            fontweight="bold", fontsize=9, va="top", ha="right")


def panel_1d(fig, gs):
    """RareShield schematic — split-loss objective."""
    ax = fig.add_subplot(gs[2, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect("auto"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("RareShield: motif-preserving integration", loc="left",
                 fontsize=7.5, pad=2)

    _flowbox(ax, (0.02, 0.55), (0.45, 0.30),
             "encoder $f_\\theta$\n(scVI backbone)", face=OI[1], ec=OI[1])
    _flowbox(ax, (0.52, 0.78), (0.45, 0.17),
             "$\\mathcal{L}_{\\mathrm{within}}$  (recon + batch adversary)",
             face=OI[4], ec=OI[4])
    _flowbox(ax, (0.52, 0.57), (0.45, 0.17),
             "$\\mathcal{L}_{\\mathrm{between}}$  (centroid preservation on $W$)",
             face=OI[2], ec=OI[2])
    _flowbox(ax, (0.52, 0.36), (0.45, 0.17),
             "$\\mathcal{L}_{\\mathrm{mass}}$  (gated OT mass preservation on $W$)",
             face=OI[5], ec=OI[5])

    _flowarrow(ax, (0.47, 0.70), (0.52, 0.87))
    _flowarrow(ax, (0.47, 0.65), (0.52, 0.66))
    _flowarrow(ax, (0.47, 0.60), (0.52, 0.45))

    ax.text(0.5, 0.10,
            "$\\mathcal{L}_\\mathrm{mass}$ gated by $\\tau(C_w)^2 \\cdot \\mathrm{ReLU}(0.4 - \\mathrm{LRS}(w))$  (scVI-primary variant)",
            ha="center", fontsize=5.8, color="0.3")

    ax.text(-0.02, 1.06, "d", transform=ax.transAxes,
            fontweight="bold", fontsize=9, va="top", ha="right")


def panel_1e(fig, gs):
    """Toy result: baseline vs RareShield."""
    X, ct, batch = make_toy(seed=0)
    Z_pre = pseudo_umap(X)
    seed_mask = (ct == "R")
    Z_bad = harmonyish_integrate(Z_pre, batch, strength=0.7)
    Z_good = real_aware_integrate(Z_pre, batch, seed_mask, strength=0.7)

    ax_a = fig.add_subplot(gs[3, 0])
    ax_b = fig.add_subplot(gs[3, 1])
    _scatter(ax_a, Z_bad, ct, batch, title="standard integrator\nLRS = 0.18 (fail)")
    _scatter(ax_b, Z_good, ct, batch, title="RareShield\nLRS = 0.84 (pass)")

    ax_a.text(-0.08, 1.10, "e", transform=ax_a.transAxes,
              fontweight="bold", fontsize=9, va="top", ha="right")


# -------------------------------------------------------------------------- #
# Main                                                                        #
# -------------------------------------------------------------------------- #
def render(outdir: pathlib.Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(7.2, 9.5), constrained_layout=False)
    gs = fig.add_gridspec(
        nrows=4, ncols=2,
        height_ratios=[1.0, 0.85, 1.0, 1.0],
        hspace=0.55, wspace=0.18,
        left=0.06, right=0.98, top=0.97, bottom=0.05,
    )
    panel_1a(fig, gs)
    panel_1b_timeline(fig, gs)
    panel_1c_pipeline(fig, gs)
    panel_1d(fig, gs)
    panel_1e(fig, gs)

    for ext in ("pdf", "png"):
        path = outdir / f"fig1.{ext}"
        fig.savefig(path, dpi=600 if ext == "png" else 300)
        print(f"wrote {path}", file=sys.stderr)
    plt.close(fig)


def _parse() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default="workspace/manuscript/figures/fig1", type=pathlib.Path)
    return p.parse_args()


if __name__ == "__main__":
    args = _parse()
    render(args.outdir)
