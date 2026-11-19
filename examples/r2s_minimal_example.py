"""Minimal R2S example (CPU, well under 1 minute after the one-time downloads).

Scores the responses of one donor's CD14 monocytes to IFN-beta, IFN-gamma and TNF-alpha (GSE181897; bundled in
examples/data) against the 90-cytokine single-cell reference, and prints the top signals.

First run: downloads the signature database r2s_signatures_v1 (~150 MB) into ~/.cache/response2signal (or set
R2S_REFERENCE to an unpacked copy) and the published CytoSig matrix (~4 MB) from its official repository.
"""
from pathlib import Path

import pandas as pd

import response2signal as r2s

resp = pd.read_csv(Path(__file__).parent / "data" / "gse181897_cd14_mono_donor_response.csv.gz", index_col=0)
print(f"response: {resp.shape[0]} genes x {resp.shape[1]} samples ({', '.join(resp.columns)})")

ref = r2s.load_reference()

# R2S (hierarchical; the default) over all 90 reference cytokines
z = r2s.score(resp, receiver="single_cell", reference=ref, output="wide")  # signals x samples
res = r2s.to_long(z)  # sample, signal, z, rank
print("\nR2S-hierarchical, top 3 of 90 signals")
print(r2s.top_signals(res, k=3).to_string(index=False))

# the same z-scores ranked within a 10-signal candidate panel (closed-set prediction)
p10 = ["IFN-gamma", "TNF-alpha", "IL-1-beta", "IL-2", "IL-4", "IL-6", "IL-10", "GM-CSF", "LAP-TGF-beta1", "IFN-beta"]
res10 = r2s.to_long(z.loc[p10])  # same as r2s.score(..., panel=p10)
print("\nR2S-hierarchical, top 3 of the 10-signal panel")
print(r2s.top_signals(res10, k=3).to_string(index=False))

# the comparator R2S-pool (one pooled signature per signal)
pool = r2s.score(resp, receiver="single_cell", method="pool", reference=ref)
print("\nR2S-pool, top 3 of 90 signals")
print(r2s.top_signals(pool, k=3).to_string(index=False))
