# Curated analysis scripts

These scripts are the curated public-data analysis snapshot used to generate the manuscript source tables and figures. They are provided for transparency and require the public datasets listed in the manuscript. Raw public data archives are not included in this repository.

Configure `PUBLIC_BIOINF_DATABASES` and project paths before re-running data acquisition or cohort-specific preprocessing scripts.

`R_compute_all_independent_meta_v2_5.R` reproduces the supplementary six-cohort
sensitivity synthesis from archived cohort-level odds ratios and confidence
intervals. It does not refit cohort models. `R_render_v2_5_submission_figures.R`
renders the submission figures from the curated source-data tables.
