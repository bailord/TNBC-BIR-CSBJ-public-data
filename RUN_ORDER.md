# Run order

1. Prepare locked gene sets.
2. Score bulk cohorts using the shared signature-scoring utility.
3. Run I-SPY2/GSE194040 treatment-context models.
4. Run GSE25066 RMA processing/audit and chemotherapy-context models.
5. Run manual external NAC cohort mapping/meta-analysis and METABRIC survival analyses.
6. Summarize GSE266919 single-cell source localization.
7. Reproduce the six-independent-cohort sensitivity synthesis with `scripts/analysis/R_compute_all_independent_meta_v2_5.R`, using only the archived cohort-level odds ratios and confidence intervals.
8. Generate manuscript figures from the archived source-data TSV files with `scripts/analysis/R_render_v2_5_submission_figures.R`.

The sensitivity script does not refit any cohort-level logistic model. Figure 2 continues to display the five-additional-cohort primary synthesis; the six-cohort result is supplementary.
