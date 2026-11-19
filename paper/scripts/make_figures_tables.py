"""Figures 2-4 and Tables 2-3 of the paper, generated from the source tables in paper/data/.

Fig. 2   single-cell benchmark            data/single_cell_benchmark.csv
Fig. 3   bulk benchmark                   data/bulk_accuracy.csv, bulk_tests.csv, bulk_held_out_tests.csv
Fig. 4   what drives a prediction         data/interpretability_genes.csv
Table 2  tables/table2_singlecell.tex     data/single_cell_benchmark.csv
Table 3  tables/table3_bulk.tex           data/bulk_accuracy.csv, bulk_tests.csv, bulk_held_out_tests.csv

Run from the repository root:  python paper/scripts/make_figures_tables.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PAPER = Path(__file__).resolve().parents[1]
DATA = PAPER / "data"
FIG = PAPER / "figures"
TAB = PAPER / "tables"

ORANGE, BLUE, GREEN, PURPLE, GREY = "#eb6834", "#2a78d6", "#199e70", "#8a5cc2", "#8c8c8c"
INK, MUTED = "#1a1a1a", "#5c5c5c"

ORDER = [  # (label in results, display label, group)
    ("R2S-hierarchical", "R2S-hierarchical", "r2s"),
    ("R2S-pool", "R2S-pool", "r2s"),
    ("R2S-atlas", "R2S-atlas", "r2s"),
    ("CytoSig (published)", "CytoSig", "sig"),
    ("CytoSig model, atlas reference", "CytoSig model, atlas reference", "sig"),
    ("SecAct (published)", "SecAct", "sig"),
    ("Signature matching (Pearson)", "Signature matching (Pearson)", "sig"),
    ("Embedding classifier, scVI zero-shot", "scVI, zero-shot", "emb"),
    ("Embedding classifier, scVI scArches", "scVI, scArches", "emb"),
    ("Embedding classifier, PCA", "PCA", "emb"),
    ("Embedding classifier, NMF", "NMF", "emb"),
    ("Geneformer native", "Geneformer", "fm"),
    ("Geneformer PCA45", "Geneformer, 45 PCs", "fm"),
    ("scGPT native", "scGPT", "fm"),
    ("scGPT PCA45", "scGPT, 45 PCs", "fm"),
]
GCOL = {"r2s": ORANGE, "sig": BLUE, "emb": GREEN, "fm": PURPLE}
PANELS = [("P10", "10 signals"), ("P37", "37 signals"), ("P90", "90 signals")]
CHANCE = {"P10": 1 / 10, "P37": 1 / 37, "P90": 1 / 90}
DS = [("gse181897", "GSE181897 (55\u201363 donors per receiver)"),
      ("gse327907", "GSE327907 (9\u201312 donors per receiver)")]
GNAME = {"r2s": "R2S variants", "sig": "Signature-based methods", "emb": "Embedding-based classifiers",
         "fm": "Foundation-model embeddings"}


def fig2():
    S = pd.read_csv(DATA / "single_cell_benchmark.csv")
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": 0.8})
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 10.0), sharey=True)
    ylab = [d for _, d, _ in ORDER]
    ypos = np.arange(len(ORDER))[::-1]
    for i, (ds, dname) in enumerate(DS):
        for j, (p, pname) in enumerate(PANELS):
            ax = axes[i, j]
            for k, (m, _, g) in enumerate(ORDER):
                r = S[(S.dataset == ds) & (S.panel == p) & (S.method == m)]
                if r.empty:
                    ax.text(0.14, ypos[k], "not covered" if m in ("CytoSig (published)", "SecAct (published)")
                            else "not run", va="center", fontsize=8.5, color=MUTED, style="italic")
                    continue
                r = r.iloc[0]
                ax.plot([r.ci_lo, r.ci_hi], [ypos[k]] * 2, color=GCOL[g], lw=1.6, alpha=0.6)
                ax.plot(r.top1_receiver_mean, ypos[k], "o", color=GCOL[g], ms=6.5 if g == "r2s" else 5,
                        mec="white", mew=0.6, zorder=3)
                if m != "R2S-hierarchical":
                    nb = int(r.n_receivers_r2s_hier_better)
                    if nb >= 3:
                        ax.text(1.045, ypos[k] - 0.12, "*", va="center", ha="center", fontsize=11, color=INK)
            ax.axvline(CHANCE[p], color=GREY, ls=":", lw=1)
            ax.set_xlim(0, 1.08)
            ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
            ax.set_yticks(ypos)
            ax.set_yticklabels(ylab, fontsize=10)
            ax.grid(axis="x", color="#e6e6e6", lw=0.6)
            ax.set_axisbelow(True)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            if i == 0:
                ax.set_title(f"{pname}", fontsize=11, weight="bold")
            if i == 1:
                ax.set_xlabel("Top-1 accuracy (mean of 5 receivers)", fontsize=9.5)
            ax.tick_params(axis="y", length=0)
            if j == 0:
                for lab, (_, _, g) in zip(ax.get_yticklabels(), ORDER):
                    lab.set_color(GCOL[g])
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", ls="", color=GCOL[g], label=GNAME[g]) for g in GNAME]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0.0, 0.035, 1, 0.965), h_pad=6.0)
    for i, (_, dname) in enumerate(DS):
        y = axes[i, 0].get_position().y1 + 0.03
        fig.text(0.012, y, ("A", "B")[i], fontsize=15, weight="bold", va="bottom")
        fig.text(0.04, y, dname, fontsize=11.5, weight="bold", va="bottom")
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"figure2_singlecell.{ext}", dpi=300)
    plt.close(fig)


def bulk_numbers():
    B = pd.read_csv(DATA / "bulk_accuracy.csv")
    T = pd.read_csv(DATA / "bulk_tests.csv")
    H = pd.read_csv(DATA / "bulk_held_out_tests.csv")
    H = H[(H.set == "held_out_1") & (H.comparison == "R2S-hierarchical vs R2S-pool")].assign(H="H1")
    return B, T, H


def h3(st):
    """R2S-hierarchical vs published CytoSig on the 30 signals CytoSig covers (held-out sets)."""
    t = pd.read_csv(DATA / "bulk_held_out_tests.csv")
    r = t[(t.set == st) & (t.comparison.str.contains("CytoSig"))].iloc[0]
    return r.rename({"top1_r2s_hierarchical": "top1_method"})


def fig3():
    B, T, H = bulk_numbers()
    sets = [("reference_loso", "Reference studies\n(LOSO, 378 studies)"), ("held_out_1", "Held-out set 1\n(121 studies)"),
            ("held_out_2", "Held-out set 2\n(135 studies)"), ("held_out_3", "Held-out set 3\n(107 studies)")]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.8), gridspec_kw={"width_ratios": [1, 1]})
    w = 0.36
    ax = axes[0]
    for i, (s, _) in enumerate(sets):
        for k, (m, col) in enumerate([("R2S-hierarchical", ORANGE), ("R2S-pool", GREY)]):
            r = B[(B.set == s) & (B.method == m) & (B.panel == "full")].iloc[0]
            x = i + (k - 0.5) * w
            ax.bar(x, r.top1, w * 0.92, color=col, label=m if i == 0 else None)
            ax.plot([x, x], [r.ci_lo, r.ci_hi], color=INK, lw=1)
    ax.set_title("Same sources, two ways of combining (41 signals)", fontsize=10.5)
    for i, (s, _) in enumerate(sets):
        t = T[(T.set == s) & (T.method == "R2S-hierarchical") & (T.comparator == "R2S-pool")].iloc[0]
        pv = H[H.H == "H1"].iloc[0].p_holm if s == "held_out_1" else t.p
        top = max(B[(B.set == s) & (B.method == m) & (B.panel == "full")].iloc[0].ci_hi
                  for m in ("R2S-hierarchical", "R2S-pool"))
        ax.text(i, top + 0.02, f"{t['diff']:+.3f}" + ("*" if pv < 0.05 else ""), ha="center", fontsize=10, color=INK)
    ax = axes[1]
    for i, (s, _) in enumerate(sets):
        if s == "reference_loso":
            h = B[(B.set == s) & (B.method == "R2S-hierarchical") & (B.panel == "full")].iloc[0].top1
            c = B[(B.set == s) & (B.method == "CytoSig procedure, rebuilt per fold")].iloc[0].top1
        else:
            r = h3(s)
            h, c = r.top1_method, r.top1_comparator
        if s == "reference_loso":
            tt = T[(T.set == s) & (T.method == "R2S-hierarchical") & (T.comparator.str.startswith("CytoSig"))].iloc[0]
            dd, pv = tt["diff"], tt.p
        else:
            dd, pv = r["diff"], (r.p_holm if s == "held_out_1" else r.p)
        ax.text(i, max(h, c) + 0.02, f"{dd:+.3f}" + ("*" if pv < 0.05 else ""), ha="center", fontsize=10, color=INK)
        ax.bar(i - 0.5 * w, h, w * 0.92, color=ORANGE, label="R2S-hierarchical" if i == 0 else None)
        ax.bar(i + 0.5 * w, c, w * 0.92, color=BLUE, hatch="//" if s == "reference_loso" else None, edgecolor="white",
               label=("CytoSig procedure rebuilt per fold (41 signals)" if s == "reference_loso" else
                      ("published CytoSig (30 signals)" if s == "held_out_1" else None)))
    ax.set_title("R2S-hierarchical vs CytoSig", fontsize=10.5)
    for ax, lab in zip(axes, "AB"):
        ax.set_xticks(range(len(sets)))
        labs = [n for _, n in sets] if lab == "A" else \
            [n if s == "reference_loso" else n.split("\n")[0] + f"\n({int(h3(s).n_studies)} studies)" for s, n in sets]
        ax.set_xticklabels(labs, fontsize=10)
        ax.set_ylim(0, 0.95)
        ax.set_ylabel("Top-1 accuracy (study mean)", fontsize=9.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.legend(frameon=False, fontsize=9.5, loc="upper left", ncol=1)
        ax.text(-0.1, 1.06, lab, transform=ax.transAxes, fontsize=14, weight="bold")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"figure3_bulk.{ext}", dpi=300)
    plt.close(fig)


def fig4():
    G = pd.read_csv(DATA / "interpretability_genes.csv")
    cases = [("GSE181897", "cd14_mono", "IFN-beta", "IFN-$\\beta$\nCD14 monocytes, GSE181897"),
             ("GSE181897", "cd14_mono", "IFN-gamma", "IFN-$\\gamma$\nCD14 monocytes, GSE181897"),
             ("GSE181897", "cd14_mono", "TNF-alpha", "TNF-$\\alpha$\nCD14 monocytes, GSE181897"),
             ("held_out_1 bulk", "all", "TNF", "TNF-$\\alpha$\nbulk receivers, held-out set 1")]
    fig, axes = plt.subplots(1, 4, figsize=(13.5, 4.4))
    for ax, (case, ct, sig, t) in zip(axes, cases):
        g = G[(G.case == case) & (G.cell_type == ct) & (G.signal == sig)].sort_values("rank").head(10)
        v = g.contribution.to_numpy() / g.contribution.max()
        cols = [ORANGE if x else "#f5c4a8" for x in g.in_signature_top20]
        ax.barh(np.arange(len(g))[::-1], v, color=cols)
        ax.set_yticks(np.arange(len(g))[::-1])
        ax.set_yticklabels(g.gene, fontsize=10, style="italic")
        ax.set_title(t, fontsize=10)
        ax.set_xlim(0, 1.05)
        ax.set_xlabel("Contribution to $\\hat{\\gamma}_c$ (scaled)", fontsize=10)
        ax.tick_params(axis="y", length=0)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    for ax, lab in zip(axes, "ABCD"):
        ax.text(-0.42, 1.13, lab, transform=ax.transAxes, fontsize=14, weight="bold")
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=ORANGE, label="also among the 20 largest weights of the consensus signature"),
                        Patch(color="#f5c4a8", label="not among them")],
               loc="lower center", ncol=2, frameon=False, fontsize=10)
    fig.tight_layout(rect=(0, 0.07, 1, 1), w_pad=3.0)
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"figure4_interpretability.{ext}", dpi=300)
    plt.close(fig)


def table2():
    S = pd.read_csv(DATA / "single_cell_benchmark.csv")
    L = [r"\begin{table}[!ht]", r"\centering",
         r"\caption{\textbf{Single-cell benchmark: top-1 accuracy (mean of 5 receivers) when the applied signal is "
         r"chosen among 10, 37 or 90 candidate signals.} Superscripts: number of receivers (of 5) in which "
         r"R2S-hierarchical was more accurate at Holm-adjusted $p<0.05$ (donor-paired sign-flip test); no method "
         r"was significantly more accurate than R2S-hierarchical in any receiver. Bold: highest accuracy in each column. SecAct covers 7 of the 10 signals of the smallest panel and is scored on the units and candidates it covers. --: the method's published reference does not cover the panel; n.r.: not "
         r"run (scArches would need a new query model). Per-receiver values and 95\% intervals of the receiver means: Supplementary Table S2.}",
         r"\label{tab:singlecell}", r"\small", r"\begin{tabular}{lcccccc}", r"\toprule",
         r" & \multicolumn{3}{c}{GSE181897} & \multicolumn{3}{c}{GSE327907} \\",
         r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}",
         r"Method & 10 & 37 & 90 & 10 & 37 & 90 \\", r"\midrule"]
    groups = {"r2s": "R2S variants", "sig": "Signature-based methods", "emb": "Embedding-based classifiers",
              "fm": "Foundation-model embeddings"}
    last = None
    for m, disp, g in ORDER:
        if g != last:
            if last is not None:
                L.append(r"\addlinespace")
            L.append(r"\multicolumn{7}{l}{\emph{" + groups[g] + r"}} \\")
            last = g
        cells = []
        for ds, _ in DS:
            for p, _ in PANELS:
                r = S[(S.dataset == ds) & (S.panel == p) & (S.method == m)]
                if r.empty:
                    cells.append("--" if m in ("CytoSig (published)", "SecAct (published)") else "n.r.")
                    continue
                r = r.iloc[0]
                v = f"{r.top1_receiver_mean:.3f}"
                col = S[(S.dataset == ds) & (S.panel == p) & S.method.isin([o[0] for o in ORDER])]
                if round(r.top1_receiver_mean, 3) == round(col.top1_receiver_mean.max(), 3):
                    v = r"\textbf{" + v + "}"
                if m != "R2S-hierarchical" and int(r.n_receivers_r2s_hier_better) > 0:
                    v += r"$^{" + str(int(r.n_receivers_r2s_hier_better)) + "}$"
                cells.append(v)
        L.append(f"{disp} & " + " & ".join(cells) + r" \\")
    L += [r"\addlinespace", r"Chance & 0.100 & 0.027 & 0.011 & 0.100 & 0.027 & 0.011 \\",
          r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    (TAB / "table2_singlecell.tex").write_text("\n".join(L) + "\n")


def pfmt(p):
    return "$<$0.001" if p < 0.001 else f"{p:.3f}"


def table3():
    B, T, H = bulk_numbers()
    f = lambda x: f"{x:.3f}"  # noqa: E731
    sets = [("reference_loso", "Reference studies (LOSO)", ""), ("held_out_1", "Held-out set 1", ""),
            ("held_out_2", "Held-out set 2", ""), ("held_out_3", "Held-out set 3", "")]
    L = [r"\begin{table}[!ht]", r"\centering",
         r"\caption{\textbf{Bulk receivers from many tissues (top-1 accuracy, study mean).} (a) R2S-hierarchical vs "
         r"R2S-pool on the 41 signals with at least three source studies. (b) R2S-hierarchical vs CytoSig. The "
         r"published CytoSig signatures were built from the reference studies, so there CytoSig's procedure was "
         r"rebuilt within each leave-one-study-out fold; in the held-out sets (studies outside the reference), "
         r"published CytoSig is compared on the 30 signals it covers. Study-paired sign-flip tests. The two "
         r"comparisons on held-out set 1 were pre-registered as one family (Holm). Differences were computed before "
         r"rounding. "
         r"LOSO: leave one study out. 95\% confidence intervals: Supplementary Table S4.}",
         r"\label{tab:bulk}", r"\footnotesize", r"\setlength{\tabcolsep}{4pt}", r"\begin{tabular}{lrlll}", r"\toprule",
         r"\multicolumn{5}{l}{\emph{(a) Same sources, two ways of combining them (41 signals)}} \\",
         r"Set & Studies & R2S-hierarchical & R2S-pool & Difference ($p$) \\", r"\midrule"]
    for s, name, status in sets:
        h = B[(B.set == s) & (B.method == "R2S-hierarchical") & (B.panel == "full")].iloc[0]
        p = B[(B.set == s) & (B.method == "R2S-pool") & (B.panel == "full")].iloc[0]
        t = T[(T.set == s) & (T.method == "R2S-hierarchical") & (T.comparator == "R2S-pool")].iloc[0]
        pv = f"Holm {H[H.H == 'H1'].iloc[0].p_holm:.3f}" if s == "held_out_1" else f"{t.p:.3f}"
        L.append(f"{name} & {int(h.n_studies)} & {f(h.top1)} & "
                 f"{f(p.top1)} & {t['diff']:+.3f} ({pv}) \\\\")
    L += [r"\midrule", r"\multicolumn{5}{l}{\emph{(b) R2S-hierarchical vs CytoSig}} \\",
          r"Set & Studies & R2S-hierarchical & CytoSig & Difference ($p$) \\", r"\midrule"]
    for s, name, status in sets:
        if s == "reference_loso":
            h = B[(B.set == s) & (B.method == "R2S-hierarchical") & (B.panel == "full")].iloc[0]
            c = B[(B.set == s) & (B.method == "CytoSig procedure, rebuilt per fold")].iloc[0]
            t = T[(T.set == s) & (T.method == "R2S-hierarchical") & (T.comparator.str.startswith("CytoSig"))].iloc[0]
            L.append(f"{name} & {int(h.n_studies)} & {f(h.top1)} & {f(c.top1)} & "
                     f"{t['diff']:+.3f} ({pfmt(t.p)}) \\\\")
        else:
            r = h3(s)
            pv = f"Holm {r.p_holm:.3f}" if s == "held_out_1" else pfmt(r.p)
            L.append(f"{name} & {int(r.n_studies)} & {f(r.top1_method)} & "
                     f"{f(r.top1_comparator)} & {r['diff']:+.3f} ({pv}) \\\\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    (TAB / "table3_bulk.tex").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    fig2()
    fig3()
    fig4()
    table2()
    table3()
    print("done")
