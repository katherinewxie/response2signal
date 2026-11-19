"""Table 1 of the paper (datasets and their roles), written from paper/data/datasets.csv.

Run from the repository root:  python paper/scripts/make_table1.py
"""
from __future__ import annotations

import csv
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]

CAPTION = (r"\caption{\textbf{Datasets and their roles.} Each dataset was assigned its role before it was scored." "\n"
           r"The single-cell cohorts were used while the method was designed, so results on them are descriptive; held-out"
           "\n" r"set 1 was built and frozen after the method was fixed. Held-out sets contain only studies outside the "
           "reference.\n" r"LOSO: leave one study out.}")


def main():
    rows = list(csv.DictReader(open(PAPER / "data" / "datasets.csv")))
    L = [r"\begin{table}[!ht]", r"\centering", CAPTION, r"\label{tab:datasets}", r"\small",
         r"\begin{tabular}{p{2.8cm}p{2.6cm}p{2.6cm}p{3.3cm}}", r"\toprule", r"Dataset & Data & Size & Role \\",
         r"\midrule"]
    block = rows[0]["block"]
    for r in rows:
        if r["block"] != block:
            L.append(r"\midrule")
            block = r["block"]
        name = r["dataset"] + (rf"~\cite{{{r['citation_key']}}}" if r["citation_key"] else "")
        L.append(f"{name} & {r['data']} & {r['size']} & {r['role']} \\\\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    out = PAPER / "tables"
    out.mkdir(parents=True, exist_ok=True)
    (out / "table1_datasets.tex").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
