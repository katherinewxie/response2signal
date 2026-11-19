"""Response2Signal (R2S): which applied extracellular signal best explains a receiver population's response?

Quick start::

    import response2signal as r2s
    ref = r2s.load_reference()                      # signature database r2s_signatures_v1 (downloaded once)
    res = r2s.score(lfc, receiver="single_cell")    # lfc: genes x samples log-fold change vs matched control
    r2s.top_signals(res, k=3)

``r2s.score`` defaults to R2S-hierarchical (per-study signatures combined by a two-level ridge; z-score of the
shared signal effect). ``method="pool"`` gives R2S-pool (one pooled signature per signal) and ``method="atlas"``
R2S-atlas (the atlas signature of the receiver's own cell type). ``receiver="bulk"`` uses the bulk reference.
``r2s.pseudobulk_response`` builds matched-control responses from an AnnData.

"""
from .core import score, to_long, top_signals
from .reference import Reference, fetch_cytosig, load_reference
from .response import pseudobulk_response

__version__ = "2.0.0"
__all__ = ["score", "to_long", "top_signals", "load_reference", "Reference", "fetch_cytosig",
           "pseudobulk_response", "__version__"]
