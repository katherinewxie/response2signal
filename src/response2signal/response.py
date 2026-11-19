"""Receiver responses from single-cell data: matched-control pseudobulk log-fold changes."""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp


def pseudobulk_response(adata, sample_key: str, condition_key: str, control: str, layer: str | None = "counts",
                        cell_type_key: str | None = None, cell_type: str | None = None, min_cells: int = 20,
                        target_sum: float = 1e4):
    """Per-sample response of one receiver population: mean log1p(CP10k) of treated minus matched control cells.

    For each sample (donor, well or replicate; ``sample_key``) and each non-control condition, the response is
    mean_cells log1p(counts / total_counts * target_sum) of the condition's cells minus the same mean over the
    sample's own ``control`` cells. Total counts are taken over all
    genes of ``adata``, so pass the full gene matrix, not a highly-variable subset. Units with fewer than
    ``min_cells`` cells in the condition or in the control are dropped.

    adata : AnnData with raw counts in ``adata.layers[layer]`` (or ``adata.X`` when ``layer`` is None or absent).
    cell_type_key / cell_type : optionally keep only the cells of one receiver cell type.
    Returns (response, units): response = genes x units DataFrame (columns "<sample>|<condition>"); units =
    DataFrame with sample, condition, n_treated, n_control per column.
    """
    obs = adata.obs
    keep = np.ones(adata.n_obs, dtype=bool)
    if cell_type_key is not None:
        if cell_type is None:
            raise ValueError("cell_type_key needs cell_type")
        keep &= (obs[cell_type_key].astype(str) == str(cell_type)).to_numpy()
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    X = sp.csr_matrix(X) if not sp.issparse(X) else X.tocsr()
    X = X[np.flatnonzero(keep)]
    samples = obs[sample_key].astype(str).to_numpy()[keep]
    conds = obs[condition_key].astype(str).to_numpy()[keep]
    size = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
    L = (sp.diags(target_sum / np.where(size > 0, size, 1.0)) @ X.astype(np.float64)).tocsr()
    L.data = np.log1p(L.data)

    def mean(idx):
        return np.asarray(L[idx].mean(axis=0)).ravel()

    cols, rows, meta = [], [], []
    for s in sorted(set(samples)):
        ic = np.flatnonzero((samples == s) & (conds == str(control)))
        if len(ic) < min_cells:
            continue
        mc = mean(ic)
        for c in sorted(set(conds[samples == s]) - {str(control)}):
            it = np.flatnonzero((samples == s) & (conds == c))
            if len(it) < min_cells:
                continue
            rows.append(mean(it) - mc)
            cols.append(f"{s}|{c}")
            meta.append({"unit": cols[-1], "sample": s, "condition": c, "n_treated": len(it), "n_control": len(ic)})
    genes = [str(g) for g in adata.var_names]
    resp = pd.DataFrame(np.stack(rows, axis=1) if rows else np.zeros((len(genes), 0)), index=genes, columns=cols)
    return resp, pd.DataFrame(meta, columns=["unit", "sample", "condition", "n_treated", "n_control"])
