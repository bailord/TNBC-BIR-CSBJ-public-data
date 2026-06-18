# Reproducibility appendix

The analysis uses prespecified locked TNBC-BIR signatures. No cohort-specific outcome-derived gene selection, coefficient reweighting, cutoff optimization, or classifier training is introduced in this submission package.

Scores are computed within each dataset or analysis subset after gene harmonization and signature coverage checks. Cross-cohort comparisons are made at the model-output level rather than by merging expression matrices across platforms.

GSE25066_RMA is retained as the chemotherapy-context negative cohort. GSE163882 is interpreted as supportive but proliferation-sensitive. Public single-cell datasets are interpreted as source-localization evidence rather than response validation.

The submission includes source-data tables for each main figure and supplementary tables for model coefficients, cohort context, signature definitions, and public scRNA source details.
