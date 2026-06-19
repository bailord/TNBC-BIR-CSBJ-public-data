# Suggested run order

1. Review signature definitions in `source_data/supplementary_tables/`.
2. Configure local paths and public data locations as described in `scripts/analysis/README.md`.
3. Run cohort-specific acquisition/preprocessing scripts only when rebuilding from public data.
4. Run scoring/modeling scripts, then figure-generation scripts.
5. Compare regenerated source tables and figures with the deposited source-data files.
