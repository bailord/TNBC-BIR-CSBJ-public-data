#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(metafor))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: compute_all_independent_meta_v2_5.R <figure2_source.tsv> <output_dir>")
}

source_path <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source <- read.delim(source_path, check.names = FALSE, stringsAsFactors = FALSE)
cohorts <- c(
  "GSE25066_RMA",
  "GSE41998",
  "GSE20194",
  "GSE22226",
  "GSE32646",
  "GSE163882"
)

rows <- source[
  source$panel == "A_pCR_summary" &
    source$dataset %in% cohorts &
    source$score_or_axis == "APC_CXCL9_axis",
  ,
  drop = FALSE
]

if (nrow(rows) != length(cohorts) || !setequal(rows$dataset, cohorts)) {
  stop("Frozen Figure 2 source does not contain exactly one eligible row per independent cohort")
}

rows <- rows[match(cohorts, rows$dataset), , drop = FALSE]
if (any(rows$estimate <= 0 | rows$conf_low <= 0 | rows$conf_high <= 0)) {
  stop("Odds ratios and confidence limits must be positive")
}

rows$yi <- log(rows$estimate)
rows$sei <- (log(rows$conf_high) - log(rows$conf_low)) / (2 * 1.96)

fit_random <- metafor::rma.uni(
  yi = rows$yi,
  sei = rows$sei,
  method = "REML"
)
fit_fixed <- metafor::rma.uni(
  yi = rows$yi,
  sei = rows$sei,
  method = "FE"
)

result <- data.frame(
  source_context = "sensitivity synthesis including the principal independent chemotherapy-context cohort",
  meta_set = "all_independent_primary_TNBC_RMA",
  include_GSE194040 = FALSE,
  include_GSE25066 = TRUE,
  k = nrow(rows),
  fixed_estimate = unname(exp(coef(fit_fixed))),
  fixed_conf_low = unname(exp(fit_fixed$ci.lb)),
  fixed_conf_high = unname(exp(fit_fixed$ci.ub)),
  fixed_p_value = unname(fit_fixed$pval),
  random_estimate = unname(exp(coef(fit_random))),
  random_conf_low = unname(exp(fit_random$ci.lb)),
  random_conf_high = unname(exp(fit_random$ci.ub)),
  random_p_value = unname(fit_random$pval),
  tau2 = unname(fit_random$tau2),
  I2 = unname(fit_random$I2),
  Q = unname(fit_random$QE),
  Q_p_value = unname(fit_random$QEp),
  included_datasets = paste(cohorts, collapse = ";"),
  interpretation = "statistically inconclusive sensitivity synthesis; not validation",
  stringsAsFactors = FALSE
)

expected <- c(
  random_estimate = 1.145889612367016,
  random_conf_low = 0.926950164914252,
  random_conf_high = 1.416541097279049,
  random_p_value = 0.208105060718329,
  I2 = 16.551235272994663,
  Q_p_value = 0.373709085521562
)
actual <- unlist(result[1, names(expected)])
if (any(abs(actual - expected) > 1e-10)) {
  stop("Frozen six-cohort REML result did not reproduce the locked expected values")
}

input_export <- rows[, c(
  "dataset", "estimate", "conf_low", "conf_high", "yi", "sei",
  "n_total", "n_event_or_pCR", "model_name", "score_or_axis"
)]
names(input_export)[names(input_export) == "estimate"] <- "odds_ratio"
names(input_export)[names(input_export) == "conf_low"] <- "ci_low"
names(input_export)[names(input_export) == "conf_high"] <- "ci_high"
names(input_export)[names(input_export) == "yi"] <- "log_odds_ratio"
names(input_export)[names(input_export) == "sei"] <- "standard_error"

write.table(
  result,
  file = file.path(output_dir, "all_independent_pCR_meta_sensitivity.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
write.table(
  input_export,
  file = file.path(output_dir, "all_independent_pCR_meta_inputs.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

session <- capture.output(sessionInfo())
writeLines(session, file.path(output_dir, "all_independent_pCR_meta_sessionInfo.txt"))

cat(sprintf(
  "k=%d; REML OR %.12f (%.12f-%.12f); p=%.12f; I2=%.12f; Qp=%.12f\n",
  result$k,
  result$random_estimate,
  result$random_conf_low,
  result$random_conf_high,
  result$random_p_value,
  result$I2,
  result$Q_p_value
))
