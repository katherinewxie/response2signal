"""Figure 1: the R2S method.

Diagram drawn from constants; the only numbers shown are reference sizes, read from paper/data/figure1_counts.json.

Run from the repository root:  python paper/scripts/make_figure1.py

Panels:
  (A) the receiver-first task: a receiver population's response to an unknown signal, against matched control
  (B) the reference: per-study signatures from three kinds of source
  (C) the hierarchical ridge and the score
  (D) the R2S variants and when to use them
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402

PAPER = Path(__file__).resolve().parents[1]
OUT_DIR = PAPER / "figures"
COUNTS = PAPER / "data" / "figure1_counts.json"

BLUE, ORANGE, GREEN, PURPLE = "#2a78d6", "#eb6834", "#199e70", "#8a5cc2"
NEUTRAL_FC, NEUTRAL_EC = "#eff3f8", "#4a6899"
INK, MUTED = "#1a1a1a", "#5c5c5c"
TITLE_FS, BOX_FS, SMALL_FS = 11.5, 9.6, 9.2


def box(ax, x, y, w, h, text="", fc=NEUTRAL_FC, ec=NEUTRAL_EC, fs=BOX_FS, weight="normal", lw=1.3, color=INK,
        ha="center"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.015", linewidth=lw,
                                edgecolor=ec, facecolor=fc, zorder=2))
    if text:
        tx = x + w / 2 if ha == "center" else x + 0.012
        ax.text(tx, y + h / 2, text, ha=ha, va="center", fontsize=fs, color=color, weight=weight, zorder=3,
                linespacing=1.4)


def arrow(ax, x0, y0, x1, y1, color=NEUTRAL_EC, lw=1.6):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=13, linewidth=lw,
                                 color=color, zorder=4))


def title(ax, x, y, letter, text):
    ax.text(x, y, f"{letter}", fontsize=TITLE_FS + 2, weight="bold", color=INK, va="top")
    ax.text(x + 0.028, y - 0.002, text, fontsize=TITLE_FS, weight="bold", color=INK, va="top")


def heat(ax, x, y, w, h, vals, cmap="RdBu_r"):
    ax.imshow(vals, extent=(x, x + w, y, y + h), cmap=cmap, vmin=-2, vmax=2, aspect="auto", zorder=3,
              interpolation="nearest")
    ax.add_patch(Rectangle((x, y), w, h, fill=False, ec=MUTED, lw=0.8, zorder=4))


def main():
    c = json.load(open(COUNTS))
    rng = np.random.default_rng(1)
    fig = plt.figure(figsize=(13.0, 8.6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # ---------------- (A) task ----------------
    title(ax, 0.015, 0.985, "A", "Receiver-first task")
    box(ax, 0.03, 0.80, 0.12, 0.09, "Receiver cells\nexposed to an\nunknown signal", fc="#fdf0e6", ec=ORANGE)
    box(ax, 0.03, 0.68, 0.12, 0.09, "Matched control\n(same donor or\nstudy)")
    arrow(ax, 0.155, 0.80, 0.20, 0.785)
    arrow(ax, 0.155, 0.725, 0.20, 0.765)
    heat(ax, 0.205, 0.67, 0.025, 0.22, rng.normal(size=(22, 1)) * 1.2)
    ax.text(0.237, 0.78, "response $y$:\nlog fold\nchange per\ngene", ha="left", va="center",
            fontsize=SMALL_FS, color=MUTED)
    ax.text(0.03, 0.645, "Which applied signal best explains $y$?", fontsize=BOX_FS, color=INK, va="top",
            weight="bold")

    # ---------------- (B) reference ----------------
    title(ax, 0.30, 0.985, "B", "Reference: one signature per signal per study")
    srcs = [
        (f"Single-cell perturbation atlas\n(Human Cytokine Dictionary;\n{c['atlas_signals']} signals, "
         f"{c['atlas_donors']} donors)", BLUE),
        (f"Bulk perturbation compendium\nre-derived from GEO\n({c['bulk_studies']} studies, "
         f"{c['bulk_experiments']} experiments)", GREEN),
        (f"Published CytoSig\nsignatures\n({c['cytosig_signals']} signals)", PURPLE),
    ]
    for i, (t, col) in enumerate(srcs):
        y = 0.80 - i * 0.105
        box(ax, 0.31, y, 0.20, 0.09, t, fc="white", ec=col, fs=SMALL_FS + 0.4)
        arrow(ax, 0.515, y + 0.045, 0.565, 0.75, color=col, lw=1.3)
    # matrix: three signals, each with several study columns
    x0, y0, hh = 0.57, 0.62, 0.27
    groups = [("IFN-$\\beta$", [BLUE, GREEN, GREEN, GREEN, PURPLE]), ("TNF-$\\alpha$", [BLUE, GREEN, GREEN, PURPLE]),
              ("IL-4", [BLUE, GREEN, PURPLE])]
    cw, gap = 0.016, 0.014
    x = x0
    for name, cols in groups:
        base = rng.normal(size=(22, 1))
        vals = base + 0.55 * rng.normal(size=(22, len(cols)))
        heat(ax, x, y0, cw * len(cols), hh, vals)
        for j, col in enumerate(cols):
            ax.add_patch(Rectangle((x + j * cw, y0 + hh + 0.004), cw, 0.012, color=col, zorder=4))
        ax.text(x + cw * len(cols) / 2, y0 - 0.012, name, ha="center", va="top", fontsize=SMALL_FS + 0.5)
        x += cw * len(cols) + gap
    ax.text(x + 0.002, y0 + hh / 2, "...", fontsize=14, va="center", color=MUTED)
    ax.text(x0, y0 + hh + 0.03, "signatures of each signal $c$ (genes $\\times$ sources; one column $x_{cs}$ per bulk study $s$)",
            fontsize=SMALL_FS, color=MUTED)
    box(ax, 0.83, 0.66, 0.155, 0.22,
        "Studies differ in cell\ntype, dose, time and\nplatform. R2S keeps\neach study's signature;\nCytoSig pools them\n"
        "into one signature\nper signal.", fc="white", ec=NEUTRAL_EC, fs=SMALL_FS)

    # ---------------- (C) model ----------------
    title(ax, 0.015, 0.545, "C", "Hierarchical ridge and score")
    box(ax, 0.03, 0.27, 0.47, 0.23, "", fc="white", ec=NEUTRAL_EC)
    ax.text(0.045, 0.475, "Fit the response to every signature at once:", fontsize=BOX_FS, va="top")
    ax.text(0.265, 0.43, r"$\tilde{y} = \sum_c \gamma_c\, \tilde{k}_c + \sum_c \sum_{s} \delta_{cs}\, \tilde{x}_{cs} + \varepsilon$",
            fontsize=13, ha="center", va="top")
    ax.text(0.265, 0.372, r"penalty  $\lambda\sum_c \gamma_c^2 + \rho\lambda \sum_{c,s} \delta_{cs}^2$,  $\rho = 10$",
            fontsize=11.5, ha="center", va="top")
    ax.text(0.045, 0.325, "$k_c$: consensus signature of signal $c$ across sources (shared effect $\\gamma_c$)\n"
            "$x_{cs}$: signature of study $s$ (deviation $\\delta_{cs}$, penalized ten times more strongly)",
            fontsize=SMALL_FS + 0.3, va="top", color=MUTED, linespacing=1.5)
    arrow(ax, 0.505, 0.385, 0.555, 0.385)
    box(ax, 0.56, 0.305, 0.17, 0.16, "Score per signal:\npermutation $z$-score of $\\hat{\\gamma}_c$\n(1,000 gene\npermutations)",
        fc="#fdf0e6", ec=ORANGE)
    arrow(ax, 0.735, 0.385, 0.775, 0.385)
    # ranked list
    sig = [("IFN-$\\beta$", 6.1), ("IFN-$\\alpha$1", 4.3), ("IFN-$\\gamma$", 2.0), ("IL-27", 1.1), ("TNF-$\\alpha$", 0.2)]
    for i, (n, z) in enumerate(sig):
        yy = 0.455 - i * 0.034
        ax.add_patch(Rectangle((0.84, yy - 0.011), 0.022 * z, 0.022, color=ORANGE if i == 0 else "#f5c4a8", zorder=3))
        ax.text(0.835, yy, n, ha="right", va="center", fontsize=SMALL_FS + 0.4)
    ax.text(0.84, 0.28, "ranked candidate signals\n(top-ranked = prediction)", fontsize=SMALL_FS, color=MUTED, va="top")

    # ---------------- (D) variants ----------------
    title(ax, 0.015, 0.20, "D", "R2S variants")
    rows = [
        ("R2S-hierarchical (default)", "Per-study signatures, hierarchical ridge as in (C). Blood and immune receivers: atlas, bulk "
         "and CytoSig; other tissues: bulk studies only.", ORANGE),
        ("R2S-pool", "One pooled signature per signal from the same sources; one ridge regression.", NEUTRAL_EC),
        ("R2S-atlas", "Atlas signatures of the receiver's own cell type; one ridge regression.", BLUE),
    ]
    for i, (n, d, col) in enumerate(rows):
        yy = 0.14 - i * 0.045
        box(ax, 0.03, yy - 0.017, 0.20, 0.036, n, fc="white", ec=col, fs=SMALL_FS + 0.5, weight="bold")
        ax.text(0.245, yy, d, fontsize=SMALL_FS + 0.4, va="center", color=INK)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"figure1_r2s.{ext}", dpi=300)
    print("wrote", OUT_DIR / "figure1_r2s.png")


if __name__ == "__main__":
    main()
