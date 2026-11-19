# Licences of the R2S signature database (r2s_signatures_v1)

The response2signal code is MIT-licensed. The data files carry the terms of the data they were derived from:

| File / arrays | Derived from | Terms |
|---|---|---|
| `atlas.npz` (composite and own-cell-type atlas signatures, expressed-gene sets) | Parse Biosciences "10 Million Human PBMCs in a Single Experiment" (Human Cytokine Dictionary; Oesinghaus et al., bioRxiv 2025) | **CC BY-NC 4.0** (stated on https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/). The signatures are adapted (averaged log fold changes) from that dataset; credit Parse Biosciences and Oesinghaus et al.; non-commercial use only. Commercial use needs a licence from Parse Biosciences. |
| `sc_bulk.npz`, `bulk.npz`, `columns.csv` (per-study and pooled bulk signatures) | Public NCBI GEO series (accessions in `columns.csv`) | NCBI places no restrictions on the use or distribution of GEO data, but submitters may claim rights in their data (https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html). The derived signatures are released under **CC BY 4.0** (https://creativecommons.org/licenses/by/4.0/); credit Response2Signal (https://github.com/katherinewxie/response2signal) and cite the GEO accessions. |
| Published CytoSig signatures | Jiang et al., Nat Methods 2021 | **Not included.** GPL-3.0 (https://github.com/data2intelligence/CytoSig/blob/master/LICENSE). `response2signal` downloads them from the official repository at run time (pinned commit, sha256-checked). |
