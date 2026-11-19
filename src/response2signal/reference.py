"""The R2S signature database (``r2s_signatures_v1``) and the published CytoSig signatures.

The database is a directory of compressed numpy files plus CSV metadata, distributed as a release asset of
https://github.com/katherinewxie/response2signal (it is too large for the git repository). It is looked up in
this order:

1. the ``path`` argument of :func:`load_reference`;
2. the environment variable ``R2S_REFERENCE`` (a database directory);
3. ``<cache>/r2s_signatures_v1`` (``<cache>`` = ``$R2S_HOME`` or ``~/.cache/response2signal``), downloaded there
   on first use from the release URL and checked against its sha256.

Contents (see ``LICENSE_DATA.md`` and ``manifest.json`` inside the database):

* single-cell / blood configuration: composite and own-cell-type atlas signatures of the Human Cytokine Dictionary
  (90 cytokines; CC BY-NC 4.0), per-study bulk signatures (467 study x signal columns, 33 signals), pooled bulk
  signatures;
* bulk configuration: per-study bulk signatures (516 columns, 41 signals) and pooled signatures (P) from the
  CytoSig-source GEO compendium;
* ``columns.csv``: one row per study column (signal, GEO study, experiments, platform, cell source, lineage);
  ``signals.csv``: the signals of each configuration and their ligand genes / CytoSig column.

The published CytoSig matrix (GPL-3.0; Jiang et al., Nat Methods 2021) is one of the sources of the single-cell
configuration. It is not bundled; :func:`fetch_cytosig` downloads it from the official repository (pinned commit,
sha256-checked) or reads the copy shipped with an installed ``secactpy`` (identical values).
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import tarfile
import tempfile
import urllib.request
from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd

DB_VERSION = "r2s_signatures_v1"
RELEASE_URL = "https://github.com/katherinewxie/response2signal/releases/download/{version}/{name}.tar"
# sha256 of the release tar files (None = not checked)
RELEASE_SHA256 = {
    "r2s_signatures_v1": "c56dc18510efa58371128bf7d3f7931432d4580cb5a421f7e63be44e5dea0de1",
    "r2s_signatures_v1_f64": "228d5d7537c5b7e175b68cab23975b706db46b3045a96a28b41810b5299e9f02",
}
CYTOSIG_URL = ("https://raw.githubusercontent.com/data2intelligence/CytoSig/"
               "eb741454f08ead71254b148b561b19690fb44f9c/CytoSig/signature.centroid")
CYTOSIG_SHA256 = "60a3152afe429a7dbb41fa0597ab9730141ce858af16a35f7063e64bfc239bde"
SECACTPY_CYTOSIG_SHA256 = "e72365938f9482b3c211610d873c037c222591a8325eb95f602c79e29c5ea2af"
RECEIVERS = ("cd14_mono", "nk", "cd4_naive", "cd8_naive", "b_naive")


def cache_dir() -> Path:
    return Path(os.environ.get("R2S_HOME", Path.home() / ".cache" / "response2signal"))


def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while b := f.read(1 << 24):
            h.update(b)
    return h.hexdigest()


def _download(url: str, dest: Path, sha: str | None):
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=dest.parent, delete=False) as tmp:
        with urllib.request.urlopen(url, timeout=120) as r:  # noqa: S310 (fixed https URLs)
            shutil.copyfileobj(r, tmp)
    got = _sha256(tmp.name)
    if sha is not None and got != sha:
        os.unlink(tmp.name)
        raise OSError(f"sha256 mismatch for {url}: expected {sha}, got {got}")
    os.replace(tmp.name, dest)


def download_reference(precision: str = "float32", dest: str | Path | None = None) -> Path:
    """Download and unpack the signature database into ``dest`` (default: the cache directory)."""
    name = DB_VERSION if precision == "float32" else DB_VERSION + "_f64"
    root = Path(dest) if dest is not None else cache_dir()
    tar = root / f"{name}.tar"
    _download(RELEASE_URL.format(version=DB_VERSION, name=name), tar, RELEASE_SHA256.get(name))
    with tarfile.open(tar) as t:
        for m in t.getmembers():  # members are plain files under <name>/
            if not (m.isfile() and m.name.startswith(f"{name}/") and ".." not in m.name):
                raise OSError(f"unexpected member {m.name!r} in {tar}")
        t.extractall(root)  # noqa: S202 (members checked above)
    tar.unlink()
    return root / name


def _resolve(path, precision, download) -> Path:
    if path is not None:
        return Path(path)
    if os.environ.get("R2S_REFERENCE"):
        return Path(os.environ["R2S_REFERENCE"])
    name = DB_VERSION if precision == "float32" else DB_VERSION + "_f64"
    p = cache_dir() / name
    if (p / "manifest.json").exists():
        return p
    if not download:
        raise FileNotFoundError(f"signature database not found at {p}; pass path=..., set R2S_REFERENCE, "
                                "or call load_reference(download=True)")
    return download_reference(precision)


def fetch_cytosig(source: str = "auto", cache: str | Path | None = None) -> pd.DataFrame:
    """Published CytoSig signatures (genes x 43 signals), read exactly as secactpy 0.3.1 ``load_cytosig()``.

    source: "auto" (cached file, else secactpy if installed, else download), "url" (official repository, pinned
    commit; sha256-checked) or "secactpy". Both sources hold identical values (checked when the database was built).
    """
    cdir = Path(cache) if cache is not None else cache_dir()
    f = cdir / "cytosig_signature.centroid"
    if source in ("auto", "url") and f.exists() and _sha256(f) == CYTOSIG_SHA256:
        return _read_cytosig(f, compression=None)
    if source in ("auto", "secactpy"):
        try:
            import secactpy  # noqa: F401
            p = Path(secactpy.__file__).parent / "data" / "CytoSig.tsv.gz"
            if _sha256(p) == SECACTPY_CYTOSIG_SHA256:
                return _read_cytosig(p, compression="gzip")
        except (ImportError, OSError):
            if source == "secactpy":
                raise
    if source in ("auto", "url"):
        _download(CYTOSIG_URL, f, CYTOSIG_SHA256)
        return _read_cytosig(f, compression=None)
    raise ValueError(f"unknown CytoSig source {source!r}")


def _read_cytosig(path, compression):
    df = pd.read_csv(path, sep="\t", index_col=0, compression=compression)
    df.index.name = "gene"
    df = df.astype(np.float64)
    if df.isna().any().any():
        df = df.fillna(0.0)
    return df


class Reference:
    """Lazy view of a signature database directory. Arrays are float64 in memory whatever the stored precision."""

    def __init__(self, path, cytosig: pd.DataFrame | None = None):
        self.path = Path(path)
        self.manifest = json.load(open(self.path / "manifest.json"))
        self._cytosig = cytosig

    def __repr__(self):
        return f"Reference({self.manifest['name']!r}, precision={self.manifest['precision']}, path='{self.path}')"

    def _npz(self, name):
        z = np.load(self.path / name, allow_pickle=False)
        return {k: (z[k].astype(np.float64) if z[k].dtype.kind == "f" else z[k]) for k in z.files}

    @cached_property
    def atlas(self) -> dict:
        d = self._npz("atlas.npz")
        d["genes"] = [str(g) for g in d["genes"]]
        d["classes"] = [str(c) for c in d["classes"]]
        return d

    @cached_property
    def sc_bulk(self) -> dict:
        d = self._npz("sc_bulk.npz")
        d["genes"] = [str(g) for g in d["genes"]]
        d["B_classes"] = [str(c) for c in d["B_classes"]]
        d["col_id"] = [str(c) for c in d["col_id"]]
        d["col_signal"] = [c.split("|")[1] for c in d["col_id"]]
        return d

    @cached_property
    def bulk(self) -> dict:
        d = self._npz("bulk.npz")
        d["genes"] = [str(g) for g in d["genes"]]
        d["classes"] = [str(c) for c in d["classes"]]
        d["col_id"] = [str(c) for c in d["col_id"]]
        d["col_signal"] = [c.split("|")[1] for c in d["col_id"]]
        return d

    @cached_property
    def columns(self) -> pd.DataFrame:
        return pd.read_csv(self.path / "columns.csv")

    @cached_property
    def signals(self) -> pd.DataFrame:
        return pd.read_csv(self.path / "signals.csv", keep_default_na=False)

    def cytosig(self) -> pd.DataFrame:
        if self._cytosig is None:
            self._cytosig = fetch_cytosig()
        return self._cytosig

    def cytosig_map(self) -> dict:
        s = self.signals
        s = s[(s.configuration == "single_cell") & (s.cytosig_column != "")]
        return dict(zip(s.signal, s.cytosig_column))

    def expressed_genes(self, cell_type: str) -> list[str]:
        """Atlas genes detected in >= 5% of the receiver's PBS cells (the gene rule of the manuscript's cohorts)."""
        _check_ct(cell_type)
        m = self.atlas[f"expressed_{cell_type}"]
        return [g for g, k in zip(self.atlas["genes"], m) if k]

    # ---- convenience views (genes x columns DataFrames)
    def atlas_signature(self, cell_type: str | None = None) -> pd.DataFrame:
        """Own-cell-type atlas signature of ``cell_type`` (R2S-atlas) or, with None, the composite atlas signature."""
        key = "comp" if cell_type is None else f"own_{_check_ct(cell_type)}"
        return pd.DataFrame(self.atlas[key], index=self.atlas["genes"], columns=self.atlas["classes"])

    def study_signatures(self, receiver: str = "single_cell") -> pd.DataFrame:
        """Per-study signatures x_cs (raw median log-FC); column metadata in ``columns``."""
        d = self.sc_bulk if receiver == "single_cell" else self.bulk
        M = d["S"] if receiver == "single_cell" else d["X"]
        return pd.DataFrame(M, index=d["genes"], columns=d["col_id"])

    def pooled_signatures(self, receiver: str = "single_cell") -> pd.DataFrame:
        """Pooled bulk signatures: the pooled signature of the single-cell configuration or P (bulk configuration)."""
        if receiver == "single_cell":
            d = self.sc_bulk
            return pd.DataFrame(d["B_all"], index=d["genes"], columns=d["B_classes"])
        return pd.DataFrame(self.bulk["P"], index=self.bulk["genes"], columns=self.bulk["classes"])


def _check_ct(ct):
    if ct not in RECEIVERS:
        raise ValueError(f"cell_type must be one of {RECEIVERS} (atlas receivers), got {ct!r}")
    return ct


def load_reference(path: str | Path | None = None, precision: str = "float32", download: bool = True,
                   cytosig: pd.DataFrame | None = None) -> Reference:
    """Open the R2S signature database.

    path: database directory (default: ``$R2S_REFERENCE`` or the cache, downloaded on first use).
    precision: "float32" (compact release) or "float64" (exact copy used to reproduce the manuscript bit for bit).
    cytosig: optional published CytoSig DataFrame (default: :func:`fetch_cytosig` when first needed).
    """
    if precision not in ("float32", "float64"):
        raise ValueError("precision must be 'float32' or 'float64'")
    return Reference(_resolve(path, precision, download), cytosig=cytosig)
