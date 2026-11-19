"""Numerical core of R2S: ridge regression of a response on signature columns, scored by gene permutation.

Two estimators are implemented, each reproducing the code that produced the manuscript's results:

* ``secact_ridge_z``: the CytoSig/SecAct ridge as implemented in secactpy 0.3.1
  (``secact_activity_inference(..., is_differential=True, is_group_sig=False)``, numpy backend). Used by
  R2S-pool and R2S-atlas. It is re-implemented here with numpy/scipy/pandas only; secactpy is not needed.
* ``hier_ridge_z``: the two-level ("hierarchical") ridge of R2S (manuscript Methods). Columns = [shared signal columns k_c, per-study columns x_cs]; penalties lambda_gamma on the
  shared coefficients gamma_c and rho * lambda_gamma on the study deviations delta_cs; score = permutation
  z-score of gamma_c.

Permutation generators (the ``rng`` argument):

* ``"srand"``     C-library ``srand(seed)`` / ``rand()`` with the Fisher-Yates shuffle of R SecAct, permutations
                  cumulative (each shuffles the previous one). This is secactpy's default. ``rand()`` is emulated
                  here exactly as GNU libc's TYPE_3 additive-feedback generator, so the result is the same on
                  every platform (secactpy itself calls the platform's libc, which differs on macOS).
* ``"legacy"``    ``numpy.random.RandomState(seed).permutation`` (secactpy ``rng_method="numpy"``).
* ``"generator"`` ``numpy.random.default_rng(seed).permutation`` (the hierarchical ridge).
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy import linalg

EPS = 1e-12  # secactpy ridge.EPS
RAND_MAX = 2147483647  # glibc RAND_MAX
_M32 = 0xFFFFFFFF
RNGS = ("srand", "legacy", "generator")


# ----------------------------------------------------------------------------------------------- glibc rand()
def _glibc_state(seed: int) -> np.ndarray:
    """The 31-word state window of glibc random_r (TYPE_3) after srand(seed) and the 310 discarded draws."""
    seed &= _M32
    if seed == 0:  # glibc: "We must make sure the seed is not 0."
        seed = 1
    w = seed - (1 << 32) if seed >= (1 << 31) else seed
    r = [0] * 344
    r[0] = w
    for i in range(1, 31):
        hi = int(w / 127773)  # C integer division truncates toward zero
        lo = w - hi * 127773
        w = 16807 * lo - 2836 * hi
        if w < 0:
            w += 2147483647
        r[i] = w
    for i in range(31, 34):
        r[i] = r[i - 31]
    for i in range(34, 344):
        r[i] = (r[i - 31] + r[i - 3]) & _M32
    return np.array([x & _M32 for x in r[313:344]], dtype=np.uint64)


_BLOCK_CACHE: dict[int, np.ndarray] = {}


def _block_matrix(B: int) -> np.ndarray:
    """Coefficients (mod 2^32) of the next B outputs as linear functions of the 31-word window."""
    if B not in _BLOCK_CACHE:
        coef = np.zeros((31 + B, 31), dtype=np.uint64)
        coef[:31] = np.eye(31, dtype=np.uint64)
        m = np.uint64(_M32)
        for b in range(B):
            coef[31 + b] = (coef[b] + coef[28 + b]) & m
        _BLOCK_CACHE[B] = coef[31:].copy()
    return _BLOCK_CACHE[B]


def glibc_rand(seed: int, n: int, block: int = 3100) -> np.ndarray:
    """The first ``n`` values of glibc ``rand()`` after ``srand(seed)`` (int64 array).

    r_i = r_{i-31} + r_{i-3} (mod 2^32), output r_i >> 1; generated in blocks with the linear recurrence written
    as a matrix (uint64 arithmetic wraps mod 2^64, so the low 32 bits are exact).
    """
    if n <= 0:
        return np.zeros(0, dtype=np.int64)
    w = _glibc_state(seed)
    M = _block_matrix(block)
    m = np.uint64(_M32)
    out = np.empty(((n + block - 1) // block) * block, dtype=np.uint64)
    for k in range(0, len(out), block):
        blk = (M @ w) & m
        out[k:k + block] = blk
        w = blk[-31:]
    return (out[:n] >> np.uint64(1)).astype(np.int64)


def srand_permutations(n: int, n_perm: int, seed: int = 0) -> np.ndarray:
    """Forward permutation table (n_perm x n) of R SecAct / secactpy ``CStdlibRNG``.

    Fisher-Yates: for i in 0..n-2, j = i + rand() / (RAND_MAX / (n - i) + 1), swap(a[i], a[j]); each permutation
    continues from the previous one. Row p is the array after the (p+1)-th shuffle, so that the permuted response
    is ``Y[row]``.
    """
    out = np.empty((n_perm, n), dtype=np.intp)
    if n < 2:
        out[:] = np.arange(n)
        return out
    r = glibc_rand(seed, n_perm * (n - 1)).reshape(n_perm, n - 1)
    i = np.arange(n - 1, dtype=np.int64)
    J = (i + r // (RAND_MAX // (n - i) + 1)).tolist()
    a = list(range(n))
    for p in range(n_perm):
        jp = J[p]
        for ii in range(n - 1):
            j = jp[ii]
            a[ii], a[j] = a[j], a[ii]
        out[p] = a
    return out


def legacy_permutations(n: int, n_perm: int, seed: int = 0) -> np.ndarray:
    """secactpy ``generate_inverse_permutation_table_fast`` (np.random.seed + np.random.permutation), forward form."""
    rs = np.random.RandomState(seed)
    return np.stack([rs.permutation(n) for _ in range(n_perm)]) if n_perm else np.zeros((0, n), dtype=np.intp)


def permutation_table(n: int, n_perm: int, seed: int, rng: str) -> np.ndarray:
    if rng == "srand":
        return srand_permutations(n, n_perm, seed)
    if rng == "legacy":
        return legacy_permutations(n, n_perm, seed)
    if rng == "generator":
        g = np.random.default_rng(seed)
        return np.stack([g.permutation(n) for _ in range(n_perm)])
    raise ValueError(f"rng must be one of {RNGS}, got {rng!r}")


# ----------------------------------------------------------------------------------------------- column scalings
def zs(M: np.ndarray) -> np.ndarray:
    """Column z-score over rows (ddof 1); constant columns -> 0."""
    mu = M.mean(axis=0)
    sd = M.std(axis=0, ddof=1)
    with np.errstate(all="ignore"):
        out = (M - mu) / sd
    out[:, ~(sd > 0)] = 0.0
    return out


def zcol(M: np.ndarray) -> np.ndarray:
    """z-score each column over its non-missing rows (ddof 0); NaN stays NaN; constant columns are centred only."""
    M = M.astype(np.float64).copy()
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-missing columns stay NaN
        mu = np.nanmean(M, axis=0)
        sd = np.nanstd(M, axis=0)
    sd[~(sd > 0)] = 1.0
    return (M - mu) / sd


def combine(blocks) -> np.ndarray:
    """Gene-wise mean of the available (non-missing) z-scored sources; missing in every source -> 0."""
    st = np.stack(blocks)
    n = np.sum(np.isfinite(st), axis=0)
    with np.errstate(all="ignore"):
        v = np.nansum(st, axis=0) / n
    v[n == 0] = 0.0
    return v


# ----------------------------------------------------------------------------------------------- secact ridge
def secact_ridge_z(X: pd.DataFrame, Y: pd.DataFrame, lam: float = 1e4, n_perm: int = 1000, seed: int = 0,
                   rng: str = "srand", perm: np.ndarray | None = None) -> pd.DataFrame:
    """z-scores (signatures x samples) of secactpy 0.3.1 ``secact_activity_inference(Y, sig_matrix=X,
    is_differential=True, is_group_sig=False, lambda_=lam, n_rand=n_perm, seed=seed)`` (numpy backend).

    Genes = ``Y.index`` intersected with ``X.index`` (pandas order rules, as secactpy). X and Y columns are z-scored
    (ddof 1, constant -> 0); T = (X'X + lam I)^-1 X'; beta = T Y; z = (beta - mean_perm) / sd_perm with the population
    variance of the permuted betas (secactpy). ``perm``: a precomputed forward permutation table for these genes.
    """
    common = Y.index.intersection(X.index)
    Xa = X.loc[common].astype(np.float64)
    Ya = Y.loc[common].astype(np.float64)
    Xs = ((Xa - Xa.mean()) / Xa.std(ddof=1)).fillna(0)
    Ys = ((Ya - Ya.mean()) / Ya.std(ddof=1)).fillna(0)
    Xv, Yv = Xs.values, Ys.values
    n, p = Xv.shape
    m = Yv.shape[1]
    XtX = Xv.T @ Xv
    XtX_reg = XtX + lam * np.eye(p, dtype=np.float64)
    try:
        L = linalg.cholesky(XtX_reg, lower=True)
        inv = linalg.cho_solve((L, True), np.eye(p, dtype=np.float64))
    except linalg.LinAlgError:
        inv = linalg.pinv(XtX_reg)
    T = inv @ Xv.T
    beta = T @ Yv
    if perm is None:
        perm = permutation_table(n, n_perm, seed, rng)
    aver = np.zeros((p, m))
    aver_sq = np.zeros((p, m))
    T = np.ascontiguousarray(T)
    Yv = np.ascontiguousarray(Yv)
    if m < p:  # secactpy "Y-row" branch
        for i in range(n_perm):
            bp = T @ Yv[perm[i], :]
            aver += bp
            aver_sq += bp ** 2
    else:  # secactpy "T-column" branch: T[:, inv_perm] @ Y
        rows = np.arange(n)
        for i in range(n_perm):
            invp = np.empty(n, dtype=np.intp)
            invp[perm[i]] = rows
            bp = T[:, invp] @ Yv
            aver += bp
            aver_sq += bp ** 2
    mean = aver / n_perm
    var = (aver_sq / n_perm) - (mean ** 2)
    se = np.sqrt(np.maximum(var, 0.0))
    z = np.where(se > EPS, (beta - mean) / se, 0.0)
    return pd.DataFrame(z, index=list(Xs.columns), columns=list(Y.columns))


# ----------------------------------------------------------------------------------------------- hierarchical ridge
def hier_ridge_z(K: np.ndarray, S0: np.ndarray | None, Y: np.ndarray, rho: float | None, lam: float = 1e4,
                 n_perm: int = 1000, perms: np.ndarray | None = None, seed: int = 0) -> np.ndarray:
    """Two-level ridge; returns C x U permutation z-scores of the shared coefficients gamma.

    K: genes x C shared columns (no NaN); S0: genes x nS study columns (missing = 0) or None (rho -> infinity: ridge
    on K only); Y: genes x U responses (no NaN). All columns and every response are z-scored over genes (ddof 1).
    Penalty: lam on gamma, rho * lam on the study deviations delta. ``perms``: forward permutation table
    (n_perm x genes), default ``default_rng(seed)``. z = (obs - mean_perm) / sd_perm (ddof 1); sd 0 -> z 0.
    """
    n, C = K.shape
    if S0 is None or rho is None:
        Z = zs(K)
        lamv = np.full(C, lam)
    else:
        Z = zs(np.concatenate([K, S0], axis=1))
        lamv = np.concatenate([np.full(C, lam), np.full(S0.shape[1], lam * rho)])
    M = np.linalg.solve(Z.T @ Z + np.diag(lamv), Z.T)[:C]
    Yz = zs(Y)
    obs = M @ Yz
    if perms is None:
        perms = permutation_table(n, n_perm, seed, "generator")
    s1 = np.zeros_like(obs)
    s2 = np.zeros_like(obs)
    for k in range(n_perm):
        t = M @ Yz[perms[k]]
        s1 += t
        s2 += t * t
    mean = s1 / n_perm
    var = (s2 - n_perm * mean ** 2) / (n_perm - 1)
    sd = np.sqrt(np.maximum(var, 0))
    with np.errstate(all="ignore"):
        z = (obs - mean) / sd
    return np.where(sd > 0, z, 0.0)
