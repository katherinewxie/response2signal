"""R2S scoring: which applied signal best explains a receiver's response?

``score(response, ...)`` takes one or more receiver responses (log-fold changes of treated vs matched control,
genes x samples) and returns, for every signal of the reference panel, a z-score and its rank per sample.

Methods (``method=``):

* ``"hierarchical"`` (default; **R2S**). Per-study signatures x_cs are kept separate and fitted jointly with one
  shared column k_c per signal:  y = sum_c gamma_c k_c + sum_{c,s} beta_cs x_cs + e,  beta_cs = delta_cs,
  penalised by lambda_gamma ||gamma||^2 + lambda_delta ||delta||^2 with lambda_gamma = 1e4 and
  lambda_delta = 10 lambda_gamma. The score of signal c is the permutation z-score of gamma_c (1,000 gene
  permutations). Shared column k_c:
    - single-cell / blood receivers (``receiver="single_cell"``): gene-wise mean of the
      z-scored composite atlas signature, the z-scored mean of the z-scored per-study bulk columns, and the z-scored
      published CytoSig signature (where available);
    - bulk receivers (``receiver="bulk"``): mean of the signal's z-scored per-study columns.
* ``"pool"`` (**R2S-pool**, the comparator): one pooled signature per signal from the same sources, scored with
  the CytoSig ridge (lambda 1e4, permutation z). Single-cell: the consensus k_c above with the pooled bulk
  signature in place of the study columns. Bulk: the pooled bulk signature P.
* ``"atlas"`` (**R2S-atlas**; single-cell only): the atlas signature of the receiver's own cell type (``cell_type`` required) with the CytoSig ridge.

Estimand and scope: the score ranks the *applied* perturbations of a closed reference panel by how well their
total response (direct effects plus secondary mediators in the source experiments) matches the query response.
It needs a matched control; it does not identify the signal acting directly on the receiver, and it cannot
report a signal outside the panel.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from . import _ridge
from .reference import Reference, load_reference

LAMBDA = 1e4
RHO = 10.0
METHODS = ("hierarchical", "pool", "atlas")
RECEIVERS = ("single_cell", "bulk")
DEFAULT_RNG = {"hierarchical": "generator", "pool": "srand", "atlas": "srand"}


# ----------------------------------------------------------------------------------------------- inputs
def _as_frame(response, genes) -> pd.DataFrame:
    if isinstance(response, pd.Series):
        Y = response.to_frame(name=response.name if response.name is not None else "sample")
    elif isinstance(response, pd.DataFrame):
        Y = response
    else:
        a = np.asarray(response, dtype=np.float64)
        if a.ndim == 1:
            a = a[:, None]
        if genes is None:
            raise ValueError("genes is required when response is an array")
        Y = pd.DataFrame(a, columns=[f"sample{i}" for i in range(a.shape[1])])
    if genes is not None:
        if len(genes) != Y.shape[0]:
            raise ValueError(f"len(genes) = {len(genes)} but the response has {Y.shape[0]} rows")
        Y = Y.set_axis([str(g) for g in genes], axis=0)
    Y = Y.astype(np.float64)
    Y.index = Y.index.astype(str)
    Y.columns = [str(c) for c in Y.columns]
    if Y.index.has_duplicates:
        raise ValueError("duplicate gene symbols in the response; aggregate them first")
    if Y.columns.has_duplicates:
        raise ValueError("duplicate sample names in the response")
    return Y


def _patterns(F: np.ndarray):
    """Group samples by their finite-gene pattern, in order of first appearance: [(row mask, [sample idx])]."""
    pats = {}
    for u in range(F.shape[1]):
        pats.setdefault(F[:, u].tobytes(), (F[:, u], []))[1].append(u)
    return list(pats.values())


# ----------------------------------------------------------------------------------------------- single-cell
def _sc_sources(ref: Reference, G: list[str], classes: list[str]):
    """z(A), z(C), z-scored study columns (with class index) and z(B_z) on the genes G."""
    at = ref.atlas
    gi = {g: i for i, g in enumerate(at["genes"])}
    A = at["comp"][[gi[g] for g in G]]
    Az = _ridge.zcol(A)
    C = ref.cytosig()
    cmap = ref.cytosig_map()
    Cz = np.full((len(G), len(classes)), np.nan)
    cg = {g: i for i, g in enumerate(C.index)}
    ci = np.array([cg.get(g, -1) for g in G])
    for j, c in enumerate(classes):
        if c in cmap:
            v = C[cmap[c]].to_numpy()
            Cz[ci >= 0, j] = v[ci[ci >= 0]]
    Cz = _ridge.zcol(Cz)
    sb = ref.sc_bulk
    bgi = {g: i for i, g in enumerate(sb["genes"])}
    rows = np.array([bgi.get(g, -1) for g in G])
    cpos = {c: j for j, c in enumerate(classes)}
    keep = np.array([c in cpos for c in sb["col_signal"]])
    Sx = np.full((len(G), int(keep.sum())), np.nan)
    Sx[rows >= 0] = sb["S"][rows[rows >= 0]][:, keep]
    cls = np.array([cpos[c] for c, k in zip(sb["col_signal"], keep) if k])
    nonempty = np.isfinite(Sx).any(axis=0)
    Sx, cls = Sx[:, nonempty], cls[nonempty]
    Sz = _ridge.zcol(Sx)
    Bz = np.full((len(G), len(classes)), np.nan)
    for j in range(len(classes)):
        k = cls == j
        if k.any():
            with np.errstate(all="ignore"), warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)  # genes missing in every study stay NaN
                Bz[:, j] = np.nanmean(Sz[:, k], axis=1)
    return Az, Cz, Sz, cls, _ridge.zcol(Bz)


def _sc_pool_columns(ref: Reference, G, classes, Az, Cz):
    """R2S-pool consensus: mean of z(A), z(B), z(C)."""
    sb = ref.sc_bulk
    bgi = {g: i for i, g in enumerate(sb["genes"])}
    bc = {c: i for i, c in enumerate(sb["B_classes"])}
    gi = np.array([bgi.get(g, -1) for g in G])
    B = np.full((len(G), len(classes)), np.nan)
    for j, c in enumerate(classes):
        if c in bc:
            col = sb["B_all"][:, bc[c]]
            B[gi >= 0, j] = col[gi[gi >= 0]]
    return _ridge.combine([Az, _ridge.zcol(B), Cz])


def _secact_full(Y: pd.DataFrame, S: np.ndarray, classes, n_perm, seed, rng):
    """CytoSig ridge on all non-zero signature columns; all-zero columns score -inf."""
    keep = [j for j in range(S.shape[1]) if np.any(S[:, j] != 0)]
    X = pd.DataFrame(S[:, keep], index=Y.index, columns=[classes[j] for j in keep])
    Z = _ridge.secact_ridge_z(X, Y, LAMBDA, n_perm, seed, rng)
    out = np.full((len(classes), Y.shape[1]), -np.inf)
    for j in keep:
        out[j] = Z.loc[classes[j], Y.columns].to_numpy()
    return out


def _score_single_cell(ref, Y, method, cell_type, n_perm, seed, rng, expressed_only):
    at = ref.atlas
    classes = list(at["classes"])
    allowed = set(ref.expressed_genes(cell_type)) if expressed_only else set(at["genes"])
    G = [g for g in Y.index if g in allowed]  # the response's gene order
    if len(G) < 100:
        raise ValueError(f"only {len(G)} response genes are in the reference; check gene symbols (HGNC, human)")
    Yg = Y.loc[G]
    F = np.isfinite(Yg.to_numpy())
    out = np.full((len(classes), Y.shape[1]), np.nan)
    if method == "atlas":
        gi = {g: i for i, g in enumerate(at["genes"])}
        S = at[f"own_{cell_type}"][[gi[g] for g in G]]
    else:
        Az, Cz, Sz, cls, Bzz = _sc_sources(ref, G, classes)
        if method == "hierarchical":
            K = _ridge.combine([Az, Bzz, Cz])
            S0 = np.nan_to_num(Sz, nan=0.0)
        else:
            S = _sc_pool_columns(ref, G, classes, Az, Cz)
    gen = np.random.default_rng(seed)
    for mask, us in _patterns(F):
        r = np.flatnonzero(mask)
        Yp = Yg.iloc[r, us]
        if method == "hierarchical":
            perms = _ridge.permutation_table(len(r), n_perm, seed, rng) if rng != "generator" else \
                np.stack([gen.permutation(len(r)) for _ in range(n_perm)])
            out[:, us] = _ridge.hier_ridge_z(K[r], S0[r], Yp.to_numpy(), RHO, LAMBDA, n_perm, perms)
        else:
            out[:, us] = _secact_full(Yp, S[r], classes, n_perm, seed, rng)
    return classes, out, {"n_genes": len(G)}


# ----------------------------------------------------------------------------------------------- bulk
def _score_bulk(ref, Y, method, n_perm, seed, rng):
    b = ref.bulk
    classes = list(b["classes"])
    genes = b["genes"]  # the reference's gene order
    yi = {g: i for i, g in enumerate(Y.index)}
    present = np.array([g in yi for g in genes])
    if present.sum() < 100:
        raise ValueError(f"only {int(present.sum())} response genes are in the bulk reference; check gene symbols")
    Yr = np.full((len(genes), Y.shape[1]), np.nan)
    Yr[present] = Y.to_numpy()[[yi[g] for g, k in zip(genes, present) if k]]
    F = np.isfinite(Yr)
    out = np.full((len(classes), Y.shape[1]), np.nan)
    if method == "hierarchical":
        cpos = {c: j for j, c in enumerate(classes)}
        cls = np.array([cpos[c] for c in b["col_signal"]])
        Sz = _ridge.zcol(b["X"])  # over the reference genes R
        for j in range(len(classes)):
            k = cls == j
            bad = np.isfinite(Sz[:, k]).sum(axis=1) < 3
            sub = Sz[:, k]
            sub[bad] = np.nan
            Sz[:, k] = sub
    gen = np.random.default_rng(seed)
    gnames = np.array(genes)
    for mask, us in _patterns(F):
        r = np.flatnonzero(mask)
        if len(r) < 1000:
            continue  # unscorable (< 1,000 measured reference genes)
        if method == "hierarchical":
            S0 = np.where(np.isfinite(Sz[r]), Sz[r], 0.0)
            K = np.stack([S0[:, cls == j].mean(axis=1) for j in range(len(classes))], axis=1)
            perms = _ridge.permutation_table(len(r), n_perm, seed, rng) if rng != "generator" else \
                np.stack([gen.permutation(len(r)) for _ in range(n_perm)])
            out[:, us] = _ridge.hier_ridge_z(K, S0, Yr[np.ix_(r, us)], RHO, LAMBDA, n_perm, perms)
        else:
            Yp = pd.DataFrame(Yr[np.ix_(r, us)], index=gnames[r], columns=[Y.columns[u] for u in us])
            X = pd.DataFrame(np.nan_to_num(b["P"][r], nan=0.0), index=gnames[r], columns=classes)
            Z = _ridge.secact_ridge_z(X, Yp, LAMBDA, n_perm, seed, rng)
            out[:, us] = Z.loc[classes].to_numpy()
    return classes, out, {"n_genes": int(present.sum())}


# ----------------------------------------------------------------------------------------------- public API
def score(response, genes=None, receiver: str = "single_cell", method: str = "hierarchical",
          cell_type: str | None = None, n_perm: int = 1000, seed: int = 0, rng: str | None = None,
          panel=None, expressed_only: bool = False, reference: Reference | None = None,
          output: str = "long") -> pd.DataFrame:
    """Score every reference signal for one or more receiver responses.

    Parameters
    ----------
    response : DataFrame (genes x samples), Series (one sample) or array (with ``genes``).
        Log-fold change of treated vs matched control in the same receiver population, e.g. the difference of mean
        log1p(CP10k) of treated and control cells of one donor (:func:`pseudobulk_response`), or a bulk log2 fold
        change. Human HGNC gene symbols. NaN = gene not measured in that sample.
    receiver : "single_cell" (single-cell or blood receivers; 90 atlas cytokines) or "bulk" (bulk receivers from
        any tissue; 41 signals).
    method : "hierarchical" (default; R2S), "pool" (R2S-pool) or "atlas" (R2S-atlas; single-cell only).
    cell_type : atlas receiver ("cd14_mono", "nk", "cd4_naive", "cd8_naive", "b_naive"); selects the own-cell-type
        signature of method "atlas" (required there) and the gene set of ``expressed_only``. Not used otherwise.
    expressed_only : single-cell only; keep only genes detected in >= 5% of the ``cell_type``'s atlas control (PBS)
        cells. The manuscript's cohorts used this 5% rule (on cohort-specific atlas cell subsets), so the gene sets
        are close but not identical. Default False: all response genes found in the reference.
    n_perm, seed : gene permutations and their seed (1,000 and 0 in the manuscript).
    rng : permutation generator, "generator" (numpy default_rng; default for hierarchical), "srand" (C-library
        rand as in SecAct/secactpy; default for pool and atlas) or "legacy" (numpy RandomState).
    panel : optional list of candidate signals. Every reference signal is always in the regression; ``panel`` only
        restricts the ranking (closed-set prediction among the candidates).
    output : "long" (columns sample, signal, z, rank; rank 1 = best) or "wide" (signals x samples z-scores).

    Samples with the same set of measured genes share one permutation set; sets are drawn in order of first
    appearance from one generator seeded with ``seed``. Bulk samples with < 1,000 measured reference genes are not
    scored (NaN). Signals whose single-cell signature is zero on every gene score -inf.
    """
    if receiver not in RECEIVERS:
        raise ValueError(f"receiver must be one of {RECEIVERS}")
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}")
    if method == "atlas" and (receiver != "single_cell" or cell_type is None):
        raise ValueError("method='atlas' needs receiver='single_cell' and a cell_type")
    if expressed_only and (receiver != "single_cell" or cell_type is None):
        raise ValueError("expressed_only needs receiver='single_cell' and a cell_type")
    rng = rng or DEFAULT_RNG[method]
    ref = reference if reference is not None else load_reference()
    Y = _as_frame(response, genes)
    if receiver == "single_cell":
        classes, Z, _ = _score_single_cell(ref, Y, method, cell_type, n_perm, seed, rng, expressed_only)
    else:  # bulk signatures are not cell-type-specific: cell_type is not used
        classes, Z, _ = _score_bulk(ref, Y, method, n_perm, seed, rng)
    W = pd.DataFrame(Z, index=classes, columns=Y.columns)
    W.index.name = "signal"
    if panel is not None:
        missing = [c for c in panel if c not in W.index]
        if missing:
            raise ValueError(f"signals not in the {receiver} reference: {missing}")
        W = W.loc[list(panel)]
    if output == "wide":
        return W
    return to_long(W)


def to_long(W: pd.DataFrame) -> pd.DataFrame:
    """Wide z (signals x samples) -> long table (sample, signal, z, rank); ties keep the reference order."""
    parts = []
    for s in W.columns:
        z = W[s].to_numpy()
        if not np.isfinite(z).any() and np.isnan(z).all():
            parts.append(pd.DataFrame({"sample": s, "signal": W.index, "z": z, "rank": np.nan}))
            continue
        order = np.argsort(-np.nan_to_num(z, nan=-np.inf), kind="stable")
        rank = np.empty(len(z))
        rank[order] = np.arange(1, len(z) + 1)
        parts.append(pd.DataFrame({"sample": s, "signal": W.index, "z": z, "rank": rank.astype(int)}))
    L = pd.concat(parts, ignore_index=True)
    return L.sort_values(["sample", "rank"], kind="stable").reset_index(drop=True) if len(L) else L


def top_signals(result: pd.DataFrame, k: int = 3) -> pd.DataFrame:
    """Top-k signals per sample from a long ``score`` result (columns sample, top1..topk, z1..zk)."""
    rows = []
    for s, g in result.groupby("sample", sort=False):
        g = g.dropna(subset=["rank"]).sort_values("rank").head(k)
        r = {"sample": s}
        for i, (sig, z) in enumerate(zip(g.signal, g.z), start=1):
            r[f"top{i}"], r[f"z{i}"] = sig, float(z)
        rows.append(r)
    return pd.DataFrame(rows)
