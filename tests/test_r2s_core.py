"""Unit tests of the R2S package (no data files needed; a small synthetic reference is built in tmp_path)."""
from __future__ import annotations

import ctypes
import ctypes.util
import json
import platform
import warnings

import numpy as np
import pandas as pd
import pytest

import response2signal as r2s
from response2signal import _ridge
from response2signal.reference import RECEIVERS


# ----------------------------------------------------------------------------------------------- RNG
@pytest.mark.skipif(platform.system() != "Linux" or platform.libc_ver()[0] != "glibc", reason="needs glibc")
@pytest.mark.parametrize("seed", [0, 1, 7, 2**31 + 5])
def test_glibc_rand_matches_libc(seed):
    libc = ctypes.CDLL(ctypes.util.find_library("c"))
    libc.rand.restype = ctypes.c_int
    libc.srand.argtypes = [ctypes.c_uint]
    libc.srand(ctypes.c_uint(seed))
    ref = np.array([libc.rand() for _ in range(20000)])
    assert np.array_equal(_ridge.glibc_rand(seed, 20000), ref)


def test_srand_permutations_are_cumulative_fisher_yates():
    n, k = 13, 5
    T = _ridge.srand_permutations(n, k, 0)
    r = _ridge.glibc_rand(0, k * (n - 1))
    a, pos = list(range(n)), 0
    for p in range(k):
        for i in range(n - 1):
            j = i + int(r[pos]) // (_ridge.RAND_MAX // (n - i) + 1)
            pos += 1
            a[i], a[j] = a[j], a[i]
        assert list(T[p]) == a
        assert sorted(T[p]) == list(range(n))


def test_legacy_and_generator_tables():
    rs = np.random.RandomState(3)
    assert np.array_equal(_ridge.permutation_table(50, 4, 3, "legacy"), np.stack([rs.permutation(50) for _ in range(4)]))
    g = np.random.default_rng(3)
    assert np.array_equal(_ridge.permutation_table(50, 4, 3, "generator"), np.stack([g.permutation(50) for _ in range(4)]))
    with pytest.raises(ValueError):
        _ridge.permutation_table(5, 1, 0, "other")


# ----------------------------------------------------------------------------------------------- ridge
def _toy(n=400, p=6, m=3, seed=0):
    rs = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(n)]
    X = pd.DataFrame(rs.normal(size=(n, p)), index=genes, columns=[f"s{j}" for j in range(p)])
    Y = pd.DataFrame(X.to_numpy()[:, :m] + rs.normal(size=(n, m)), index=genes, columns=[f"u{j}" for j in range(m)])
    return X, Y


@pytest.mark.parametrize("m", [2, 9])  # Y-row branch (m < p) and T-column branch (m >= p)
def test_secact_ridge_beta_and_signal(m):
    X, Y = _toy(m=min(m, 6), p=6)
    if m > 6:
        Y = pd.concat([Y, Y.add_suffix("b")], axis=1).iloc[:, :m]
    Z = _ridge.secact_ridge_z(X, Y, lam=10.0, n_perm=200, seed=0, rng="srand")
    assert Z.shape == (6, Y.shape[1]) and np.isfinite(Z.to_numpy()).all()
    for j in range(min(3, Y.shape[1])):
        assert Z.iloc[:, j].idxmax() == f"s{j % 3}"


def test_secact_branches_agree():
    X, Y = _toy(p=6, m=3)
    perm = _ridge.permutation_table(len(X), 100, 0, "srand")
    a = _ridge.secact_ridge_z(X, Y, 10.0, 100, 0, perm=perm)                     # m < p: Y-row
    Yw = pd.concat([Y] * 3, axis=1).set_axis([f"c{i}" for i in range(9)], axis=1)
    b = _ridge.secact_ridge_z(X, Yw, 10.0, 100, 0, perm=perm)                    # m >= p: T-column
    assert np.allclose(a.to_numpy(), b.to_numpy()[:, :3], atol=1e-10)


def test_hier_ridge_limits():
    rs = np.random.default_rng(1)
    K = rs.normal(size=(300, 4))
    S0 = K[:, [0, 0, 1, 2, 3]] + 0.3 * rs.normal(size=(300, 5))
    Y = K[:, [2]] + rs.normal(size=(300, 1))
    z_inf = _ridge.hier_ridge_z(K, None, Y, None, lam=1e4, n_perm=200)
    z_big = _ridge.hier_ridge_z(K, S0, Y, 1e12, lam=1e4, n_perm=200)
    assert np.allclose(z_inf, z_big, atol=1e-6)
    z = _ridge.hier_ridge_z(K, S0, Y, 10.0, lam=1e4, n_perm=200)
    assert int(np.argmax(z[:, 0])) == 2


def test_zcol_and_combine():
    M = np.array([[1.0, np.nan], [2.0, np.nan], [3.0, 5.0]])
    Z = _ridge.zcol(M)
    assert np.allclose(Z[:, 0], (M[:, 0] - 2) / np.std(M[:, 0]))
    assert np.isnan(Z[0, 1]) and Z[2, 1] == 0.0
    v = _ridge.combine([np.array([1.0, np.nan, np.nan]), np.array([3.0, 2.0, np.nan])])
    assert np.allclose(v, [2.0, 2.0, 0.0])


# ----------------------------------------------------------------------------------------------- synthetic reference
SC = ["A", "B", "C", "D", "E", "F"]
BK = ["A", "B", "C", "X"]


def _make_reference(tmp_path, n=1500, seed=0):
    rs = np.random.default_rng(seed)
    genes = np.array([f"G{i}" for i in range(n)])
    base = {c: rs.normal(size=n) for c in SC + ["X"]}
    d = tmp_path / "r2s_signatures_v1"
    d.mkdir()
    own = {ct: np.stack([base[c] + 0.3 * rs.normal(size=n) for c in SC], 1) for ct in RECEIVERS}
    np.savez_compressed(d / "atlas.npz", genes=genes, classes=np.array(SC),
                        comp=np.median(np.stack(list(own.values())), 0).astype(np.float32),
                        **{f"own_{ct}": v.astype(np.float32) for ct, v in own.items()},
                        **{f"expressed_{ct}": np.arange(n) < n - 100 for ct in RECEIVERS})
    cols, ids = [], []
    for c in SC[:4]:
        for s in range(3):
            v = base[c] + 0.8 * rs.normal(size=n)
            v[rs.random(n) < 0.05] = np.nan
            cols.append(v)
            ids.append(f"sc|{c}|GSE{s}")
    S = np.stack(cols, 1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN genes
        B = np.stack([np.nanmedian(S[:, 3 * i:3 * i + 3], 1) for i in range(4)], 1)
    np.savez_compressed(d / "sc_bulk.npz", genes=genes, B_all=B.astype(np.float32), B_classes=np.array(SC[:4]),
                        S=S.astype(np.float32), col_id=np.array(ids))
    cols, bids = [], []
    for c in BK:
        for s in range(4):
            cols.append(base[c] + 0.8 * rs.normal(size=n))
            bids.append(f"bulk|{c}|GSE{10 + s}")
    X = np.stack(cols, 1)
    P = np.stack([np.median(X[:, 4 * i:4 * i + 4], 1) for i in range(len(BK))], 1)
    np.savez_compressed(d / "bulk.npz", genes=genes, classes=np.array(BK), P=P.astype(np.float32),
                        X=X.astype(np.float32), col_id=np.array(bids))
    pd.DataFrame({"column_id": ids + bids}).to_csv(d / "columns.csv", index=False)
    pd.DataFrame({"configuration": ["single_cell"] * len(SC) + ["bulk"] * len(BK), "signal": SC + BK,
                  "cytosig_column": ["CA", "CB", "", "", "", ""] + [""] * len(BK)}).to_csv(d / "signals.csv", index=False)
    json.dump({"name": "r2s_signatures_v1", "precision": "float32"}, open(d / "manifest.json", "w"))
    cyto = pd.DataFrame({"CA": base["A"][:1000], "CB": base["B"][:1000]}, index=genes[:1000])
    return r2s.load_reference(d, cytosig=cyto), base, genes


@pytest.fixture(scope="module")
def synth(tmp_path_factory):
    return _make_reference(tmp_path_factory.mktemp("ref"))


@pytest.mark.parametrize("method,kw", [("hierarchical", {}), ("pool", {}), ("atlas", {"cell_type": "nk"}),
                                       ("hierarchical", {"cell_type": "nk", "expressed_only": True})])
def test_score_single_cell_recovers_signal(synth, method, kw):
    ref, base, genes = synth
    rs = np.random.default_rng(5)
    Y = pd.DataFrame({f"true_{c}": base[c] + 1.5 * rs.normal(size=len(genes)) for c in ["A", "C", "E"]}, index=genes)
    res = r2s.score(Y, method=method, reference=ref, n_perm=200, **kw)
    top = r2s.top_signals(res, k=1).set_index("sample").top1
    if method == "hierarchical" or method == "pool":
        assert top["true_A"] == "A" and top["true_C"] == "C"
    else:
        assert list(top) == ["A", "C", "E"]
    assert set(res.columns) == {"sample", "signal", "z", "rank"}
    assert res.groupby("sample")["rank"].apply(lambda r: sorted(r) == list(range(1, 7))).all()


@pytest.mark.parametrize("method", ["hierarchical", "pool"])
def test_score_bulk_recovers_signal_with_missing_genes(synth, method):
    ref, base, genes = synth
    rs = np.random.default_rng(6)
    Y = pd.DataFrame({f"true_{c}": base[c] + 1.5 * rs.normal(size=len(genes)) for c in BK}, index=genes)
    Y.iloc[:200, 1] = np.nan  # a second measured-gene pattern
    W = r2s.score(Y, receiver="bulk", method=method, reference=ref, n_perm=200, output="wide")
    assert list(W.idxmax()) == BK


def test_panel_restricts_ranking_only(synth):
    ref, base, genes = synth
    y = pd.Series(base["B"] + np.random.default_rng(2).normal(size=len(genes)), index=genes, name="s")
    full = r2s.score(y, reference=ref, n_perm=100, output="wide")
    sub = r2s.score(y, reference=ref, n_perm=100, output="wide", panel=["C", "B"])
    assert np.array_equal(sub.to_numpy(), full.loc[["C", "B"]].to_numpy())


def test_score_input_validation(synth):
    ref, base, genes = synth
    with pytest.raises(ValueError):
        r2s.score(np.zeros(len(genes)), reference=ref)  # array without genes
    with pytest.raises(ValueError):
        r2s.score(pd.Series(base["A"], index=genes), method="atlas", reference=ref)  # atlas needs cell_type
    with pytest.raises(ValueError):
        r2s.score(pd.Series(base["A"], index=genes), receiver="bulk", method="atlas", cell_type="nk", reference=ref)
    with pytest.raises(ValueError):
        r2s.score(pd.Series(base["A"], index=[f"x{i}" for i in range(len(genes))]), reference=ref)  # no overlap
    z = r2s.score(base["A"], genes=list(genes), reference=ref, n_perm=50, output="wide")
    assert z.shape == (6, 1)


def test_reference_views(synth):
    ref, _, genes = synth
    assert ref.atlas_signature().shape == (len(genes), 6)
    assert ref.atlas_signature("nk").shape == (len(genes), 6)
    assert ref.study_signatures("bulk").shape[1] == 16
    assert list(ref.pooled_signatures("single_cell").columns) == SC[:4]
    assert len(ref.expressed_genes("cd14_mono")) == len(genes) - 100
    with pytest.raises(ValueError):
        ref.atlas_signature("monocyte")


# ----------------------------------------------------------------------------------------------- response helper
def test_pseudobulk_response():
    ad = pytest.importorskip("anndata")
    rs = np.random.default_rng(0)
    counts = rs.poisson(2.0, size=(120, 30)).astype(np.float32)
    obs = pd.DataFrame({"donor": np.repeat(["d1", "d2"], 60),
                        "cond": np.tile(np.repeat(["PBS", "IFNG", "TNF"], 20), 2)},
                       index=[f"c{i}" for i in range(120)])
    a = ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(30)]))
    resp, units = r2s.pseudobulk_response(a, "donor", "cond", "PBS", layer=None, min_cells=20)
    assert list(resp.columns) == ["d1|IFNG", "d1|TNF", "d2|IFNG", "d2|TNF"]
    c64 = counts.astype(np.float64)
    L = np.log1p(c64 / c64.sum(1, keepdims=True) * 1e4)
    m = lambda d, c: L[((obs.donor == d) & (obs.cond == c)).to_numpy()].mean(0)  # noqa: E731
    assert np.allclose(resp["d2|TNF"], m("d2", "TNF") - m("d2", "PBS"))
    assert units.n_treated.tolist() == [20] * 4
    resp2, _ = r2s.pseudobulk_response(a, "donor", "cond", "PBS", layer=None, min_cells=21)
    assert resp2.shape[1] == 0
