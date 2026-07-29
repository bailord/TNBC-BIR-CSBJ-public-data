# Reproducibility appendix

The analysis uses literature-anchored, locked TNBC-BIR component gene sets. No
cohort-specific outcome-derived gene selection, coefficient reweighting, cutoff
optimization, or classifier training is introduced in the submitted analyses.

Gene-level values are harmonized before module scoring, and signature coverage is
checked in every cohort. Scores are computed within each cohort scoring matrix before
subtype and endpoint analyses. Cross-cohort comparisons are made at the model-output
level rather than by merging expression matrices across platforms.

The antigen-presentation summary score combines the standardized HLA class II,
CXCL9/FOLR2 macrophage, and dendritic-cell/HLA module scores. It is treated as a
correlated transcriptional summary rather than evidence that the three component
programs are statistically independent.

GSE25066_RMA is retained as the independent chemotherapy-context negative cohort.
GSE163882 is interpreted as supportive but proliferation-sensitive. GSE266919
single-cell data are used only for immune-compartment source localization.

Figure 2 reports the five-additional-NAC-cohort synthesis. A separate six-independent-
cohort sensitivity analysis adds GSE25066_RMA to those five frozen cohort-level
effects; it does not refit any cohort model. Source-data tables are provided for each
main figure, together with supplementary tables for model coefficients, cohort
context, signature definitions, coverage, and single-cell source details.
