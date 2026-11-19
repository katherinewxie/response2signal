# Response2Signal (R2S)

R2S answers a receiver-first question: **which applied extracellular signal best explains the transcriptional
response of a receiver cell population?** It compares the population's response, measured against a matched
control, with a reference of perturbation signatures. For each signal, the reference keeps one signature per source
study, drawn from a single-cell cytokine atlas, a compendium of public bulk perturbation experiments and the
published CytoSig signatures. A hierarchical ridge regression combines the study signatures.

> Katherine Xie and Jie Wu. *Response2Signal: integrating single-cell and bulk perturbation studies to identify
> extracellular signals from receiver responses.* Manuscript under review at BMC Bioinformatics.
> Department of Biological Chemistry, University of California, Irvine. Corresponding author: Jie Wu
> (jiew5@uci.edu).

![R2S overview](paper/figures/figure1_r2s.png)

*(A) The input is the response of receiver cells to an unknown signal, relative to a matched control. (B) The
reference keeps one signature per signal per study. (C) A hierarchical ridge fits all signatures at once; the score
of a signal is the permutation z-score of its shared coefficient. (D) The three R2S variants.*

## Repository contents

| Path | Contents |
|---|---|
| `src/response2signal/` | the R2S package |
| `examples/` | a minimal, runnable example with bundled data |
| `tests/` | unit tests (no data download needed) |
| `paper/` | scripts and source tables that reproduce the figures and tables of the main text (see `paper/README.md`) |

## Installation

Python 3.10 or later. The package needs only numpy, scipy and pandas.

```bash
git clone https://github.com/katherinewxie/response2signal.git
cd response2signal
pip install -e .                 # R2S
pip install -e ".[anndata]"      # + pseudobulk_response from an AnnData object
pip install -e ".[test]"         # + pytest
pytest tests                     # unit tests
```

## Quick start

```python
import pandas as pd
import response2signal as r2s

resp = pd.read_csv("examples/data/gse181897_cd14_mono_donor_response.csv.gz", index_col=0)  # genes x samples
ref = r2s.load_reference()                                   # downloads the signature database once
res = r2s.score(resp, receiver="single_cell", reference=ref)  # R2S-hierarchical (default)
print(r2s.top_signals(res, k=3))
```

The full example, which also restricts the ranking to a 10-signal candidate panel and runs R2S-pool, is

```bash
python examples/r2s_minimal_example.py
```

It scores the responses of one donor's CD14 monocytes to IFN-β, IFN-γ and TNF-α (GSE181897) and runs on a CPU in
well under a minute after the one-time downloads.

From an AnnData object of one receiver cell type with raw counts:

```python
resp, units = r2s.pseudobulk_response(adata, sample_key="donor", condition_key="condition", control="PBS")
res = r2s.score(resp, receiver="single_cell")
```

## Input and output

- **Input:** the response of a receiver population, i.e. the log fold change of treated versus **matched control**
  cells (same donor, sample or study), for one or more samples, indexed by human gene symbols.
  - Single-cell data: mean log1p(CP10k) of treated minus control cells (`r2s.pseudobulk_response`).
  - Bulk data: a log2 fold change.
- **Output:** for every signal in the reference, a z-score and its rank. The top-ranked signal is the prediction.

## Reference configurations and variants

| `receiver=` | Use for | Signals | Reference |
|---|---|---|---|
| `"single_cell"` | single-cell or blood (PBMC) receivers | 90 cytokines of the atlas | composite atlas signatures, 467 per-study bulk signatures (33 signals), published CytoSig |
| `"bulk"` | bulk receivers from any tissue | 41 signals | 516 per-study bulk signatures from the CytoSig source studies |

| `method=` | Variant | Description |
|---|---|---|
| `"hierarchical"` (default) | R2S-hierarchical | per-study signatures combined by a two-level ridge |
| `"pool"` | R2S-pool | one pooled signature per signal from the same sources |
| `"atlas"` | R2S-atlas | atlas signatures of the receiver's own cell type (single-cell only; `cell_type=` required) |

**Method in brief.** The response y and all signature columns are z-scored over genes, and

  y = Σ_c γ_c k_c + Σ_c Σ_s δ_cs x_cs + ε,  minimizing ‖ε‖² + λ Σ_c γ_c² + 10 λ Σ_{c,s} δ_cs²,  λ = 10⁴,

where k_c is the consensus signature of signal c and x_cs the signature of study s (the gene-wise median log fold
change over the study's experiments with c). The score of c is the z-score of γ_c against 1,000 gene permutations
of y. R2S-pool and R2S-atlas use one column per signal in the same ridge. The paper's Methods give the full
definitions.

## API

| Function | Purpose |
|---|---|
| `r2s.score(response, genes=None, receiver="single_cell", method="hierarchical", cell_type=None, n_perm=1000, seed=0, rng=None, panel=None, expressed_only=False, reference=None, output="long")` | z-score and rank of every signal for each sample; `output="wide"` returns a signals × samples matrix; `panel` restricts the ranking to candidate signals (the regression always uses the whole reference) |
| `r2s.top_signals(result, k=3)` | top-k signals per sample |
| `r2s.to_long(z)` | wide z-scores → long table with ranks |
| `r2s.load_reference(path=None, precision="float32")` | open the signature database |
| `r2s.fetch_cytosig()` | the published CytoSig matrix |
| `r2s.pseudobulk_response(adata, sample_key, condition_key, control)` | matched-control responses from single-cell counts |

## Signature database

The reference is the versioned database **`r2s_signatures_v1`**, distributed as a release asset of this
repository because it is too large for git:

| Asset | Size | Content |
|---|---|---|
| `r2s_signatures_v1.tar` | 148 MB | float32 signatures (default) |
| `r2s_signatures_v1_f64.tar` | 224 MB | float64 copy that reproduces the paper's z-scores exactly |

`r2s.load_reference()` downloads `r2s_signatures_v1` on first use, verifies its sha256 checksum and caches it in
`~/.cache/response2signal`. Set `R2S_HOME` to change the cache directory, or `R2S_REFERENCE` to point to an unpacked
copy. To work offline, download the asset from
`https://github.com/katherinewxie/response2signal/releases/download/r2s_signatures_v1/r2s_signatures_v1.tar`, untar
it and call `load_reference(path=...)`. The database contains the atlas signatures, the per-study and pooled bulk
signatures, and `columns.csv`, which lists the GEO study, ligand, platform and cell source of every study signature.

The published CytoSig signatures are not redistributed (GPL-3.0); `r2s.fetch_cytosig()` downloads them from the
official repository at a pinned commit and checks their sha256.

## Scope

- **R2S identifies the applied perturbation, not the directly acting signal.** A signature is the total response in
  its source experiment, including secondary effects of mediators released by other cells. A high score for signal
  c means that the response resembles the response to applying c, not that c acted directly on the receiver.
- **It is closed-set.** It ranks only the signals in the reference (`r2s.load_reference().signals`).
- **It needs a matched control.**
- **It does not infer senders,** ligand–receptor pairs or communication networks; sender-first tools (for example
  LIANA+, NicheNet or CellChat) address those questions.
- **Evaluated settings:** human PBMC receivers (CD14 monocytes, NK cells, naive CD4 and CD8 T cells, naive B cells)
  in two single-cell cohorts, and public bulk perturbation experiments from many tissues. Accuracy and its
  uncertainty are reported in the paper.

## Licences

The code is released under the MIT licence (`LICENSE`). The data have their own terms:

| Data | Terms |
|---|---|
| Atlas-derived signatures (Human Cytokine Dictionary, Parse Biosciences) | CC BY-NC 4.0, non-commercial use only |
| GEO-derived bulk signatures | CC BY 4.0; cite the GEO accessions listed in `columns.csv` |
| Published CytoSig signatures | GPL-3.0; downloaded, not redistributed |
| `examples/data/` | see `examples/data/LICENSE_DATA.txt` |

## Citation

If you use R2S, please cite the paper above (see also `CITATION.cff`).
