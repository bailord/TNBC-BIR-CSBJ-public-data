# TNBC-BIR CSBJ public-data reproducibility materials

This repository accompanies the manuscript **A prespecified computational framework reveals context-dependent APC/CXCL9/HLA-II immune-reactive states in triple-negative breast cancer**.

The repository contains curated public-data source tables, supplementary model tables, final figure files, and lightweight figure-rendering/support scripts for the CSBJ submission package. It does **not** contain raw public database downloads, protected health information, local institutional data, internal drafting packages, or manuscript-review archives.

## Contents

- `source_data/figure_source_data/`: source-data tables for the main and supplementary figures.
- `source_data/supplementary_tables/`: supplementary tables, including locked signature definitions and model summaries.
- `figures/`: final main figures, supplementary figure, and graphical abstract.
- `scripts/`: lightweight reproducibility/support scripts used for final figure packaging.
- `docs/reproducibility_appendix.md`: short reproducibility notes and cohort/source-data index.
- `repository_manifest.tsv`: file-level sizes and SHA-256 checksums.

## Data provenance

All datasets analyzed in the manuscript were obtained from public repositories or previously published resources, including GEO accessions and METABRIC-derived resources described in the manuscript and supplementary tables. Original public datasets remain governed by their source repository and publication terms.

## Analysis boundary

The TNBC-BIR framework uses prespecified locked signatures and reports context-dependent state-decomposition evidence. No outcome-derived cutoff optimization or machine-learning rescue model is included in this repository.

## License

Code and text in this repository are provided under the MIT License. Data tables derived from public repositories remain subject to the terms and licenses of the original data sources.
