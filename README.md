# TNBC-BIR CSBJ reproducibility materials

This repository archives the curated source data, locked gene-set definitions, analysis scripts, final figure source data, and supplementary tables for the CSBJ submission. It does not redistribute raw public database files, controlled data, or local hospital data.

The study maps context-dependent antigen-presentation states in triple-negative breast cancer. The antigen-presentation summary score is a correlated transcriptional summary and is not presented as a general clinical pathological-complete-response prediction assay.

## Submission release

Release `v1.6-submission` improves the column-safe line wrapping in Figure 1 and retains the corrected Figure 2 hollow-point annotation identifying the age- and stage-adjusted complete-case sensitivity analysis. No model, estimate, confidence interval, p value, sample count, or interpretation was changed. The release retains the frozen six-independent-cohort pCR sensitivity synthesis, complete figure source data, supplementary tables, and the exact scripts used to reproduce the submitted displays.

## Contents

- `scripts/analysis/`: curated cohort processing, scoring, modeling, meta-analysis, single-cell localization, and figure-generation scripts.
- `source_data/figure_source_data/`: row-level data used for Figures 1-4 and Supplementary Figure 1.
- `source_data/supplementary_tables/`: editable main tables and Supplementary Tables 1-5.
- `figures/`: final PDF, PNG, and editable SVG figures.
- `RUN_ORDER.md`: analysis and rendering sequence.

Raw public datasets should be obtained from the accessions reported in the manuscript and supplementary tables.
