#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(tidyr)
  library(yaml)
})

set.seed(42)
root <- getwd()
cfg <- yaml::read_yaml(file.path(root, "config", "project_config.yaml"))
results_dir <- file.path(root, cfg$paths$results)
figures_dir <- file.path(root, "figures")
reports_dir <- file.path(root, "reports")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", trimws(formatC(x, format = "fg", digits = digits)))
}

or_rows <- function(fit, dataset, analysis_subset, endpoint, n_total, n_event, model_name, covariates, notes) {
  sm <- coef(summary(fit))
  tibble(
    dataset = dataset,
    analysis_subset = analysis_subset,
    endpoint = endpoint,
    n_total = n_total,
    n_event_or_pCR = n_event,
    model_name = model_name,
    variable = rownames(sm),
    effect_type = "OR",
    estimate = exp(coef(fit)),
    conf_low = exp(coef(fit) - 1.96 * sm[, "Std. Error"]),
    conf_high = exp(coef(fit) + 1.96 * sm[, "Std. Error"]),
    p_value = sm[, "Pr(>|z|)"],
    covariates = covariates,
    notes = notes
  )
}

fit_model <- function(df, formula, model_name, notes) {
  fit_df <- model.frame(formula, data = df, na.action = na.omit)
  n_event <- sum(fit_df$pCR_binary == 1, na.rm = TRUE)
  fit <- glm(formula, data = df, family = binomial(), na.action = na.omit)
  or_rows(fit, "GSE163882", "primary_TNBC", "pCR", nrow(fit_df), n_event,
          model_name, deparse(formula), notes)
}

main <- function() {
  score_path <- file.path(root, "data", "processed", "external_validation_external_NAC_scores_clinical.tsv")
  qc_path <- file.path(results_dir, "external_validation_external_dataset_mapping_qc.tsv")
  coverage_path <- file.path(results_dir, "external_validation_external_signature_coverage.tsv")
  dl_path <- file.path(results_dir, "external_validation_external_download_qc.tsv")
  if (!file.exists(score_path)) stop("Missing processed Round 4 external NAC scores.")
  scores <- readr::read_tsv(score_path, show_col_types = FALSE) %>% filter(dataset == "GSE163882")
  if (!nrow(scores)) stop("GSE163882 rows are absent from processed score table.")
  qc <- readr::read_tsv(qc_path, show_col_types = FALSE) %>% filter(dataset == "GSE163882")
  coverage <- readr::read_tsv(coverage_path, show_col_types = FALSE) %>% filter(dataset == "GSE163882")
  downloads <- readr::read_tsv(dl_path, show_col_types = FALSE) %>% filter(dataset == "GSE163882")

  clinical_dist <- bind_rows(
    scores %>% count(variable = "endpoint_raw", value = endpoint_raw, name = "n"),
    scores %>% count(variable = "er_raw", value = er_raw, name = "n"),
    scores %>% count(variable = "pr_raw", value = pr_raw, name = "n"),
    scores %>% count(variable = "her2_raw", value = her2_raw, name = "n"),
    scores %>% count(variable = "primary_tnbc", value = as.character(primary_tnbc), name = "n"),
    scores %>% count(variable = "subtype_mapping_rule", value = subtype_mapping_rule, name = "n"),
    scores %>% count(variable = "cohort_note", value = cohort_note, name = "n")
  ) %>% arrange(variable, desc(n))
  readr::write_tsv(clinical_dist, file.path(results_dir, "GSE163882_clinical_variable_distribution.tsv"))

  primary <- scores %>% filter(primary_tnbc, pCR_mapped)
  expected_ok <- nrow(scores) == 222 && nrow(primary) == 90 &&
    sum(primary$pCR_binary == 1, na.rm = TRUE) == 38 &&
    sum(primary$pCR_binary == 0, na.rm = TRUE) == 52
  endpoint_ok <- all(!is.na(primary$pCR_binary)) &&
    setequal(sort(unique(primary$endpoint_raw)), sort(unique(primary$endpoint_raw[!is.na(primary$endpoint_raw)])))
  subtype_ok <- all(primary$er_negative, primary$pr_negative, primary$her2_negative, na.rm = TRUE)
  coverage_min <- min(coverage$coverage, na.rm = TRUE)
  expression_note <- unique(na.omit(scores$expression_source))[1]
  annotation_note <- unique(na.omit(scores$annotation_note))[1]

  corr <- tibble(
    dataset = "GSE163882",
    subset = "primary_TNBC",
    n = nrow(primary),
    pearson = cor(primary$APC_CXCL9_axis, primary$Proliferation, use = "pairwise.complete.obs", method = "pearson"),
    spearman = cor(primary$APC_CXCL9_axis, primary$Proliferation, use = "pairwise.complete.obs", method = "spearman"),
    notes = "Correlation calculated on primary TNBC pCR-mapped samples."
  )
  readr::write_tsv(corr, file.path(results_dir, "GSE163882_axis_proliferation_correlation.tsv"))

  models <- bind_rows(
    fit_model(primary, pCR_binary ~ scale(APC_CXCL9_axis), "axis_only",
              "Supportive if positive; audit required because this is the only statistically positive non-I-SPY2 pCR cohort."),
    fit_model(primary, pCR_binary ~ scale(APC_CXCL9_axis) + scale(Proliferation), "axis_plus_proliferation",
              "If APC axis loses significance after proliferation adjustment, interpret as supportive but not proliferation-independent."),
    fit_model(primary, pCR_binary ~ scale(Proliferation), "proliferation_only",
              "Proliferation control model."),
    fit_model(primary, pCR_binary ~ scale(TNBC_BIR_favorable_score), "aggregate_BIR_descriptive",
              "Supplementary/descriptive only; aggregate BIR is not the main claim.")
  )
  readr::write_tsv(models, file.path(results_dir, "GSE163882_axis_models_final.tsv"))

  audit <- tibble(
    dataset = "GSE163882",
    all_samples = nrow(scores),
    primary_tnbc_pcr_mapped_n = nrow(primary),
    primary_tnbc_pcr = sum(primary$pCR_binary == 1, na.rm = TRUE),
    primary_tnbc_rd = sum(primary$pCR_binary == 0, na.rm = TRUE),
    endpoint_source = "endpoint_raw from Round 4 parsed clinical table",
    endpoint_values = paste(sort(unique(primary$endpoint_raw)), collapse = ";"),
    receptor_source = "er_raw/pr_raw/her2_raw from Round 4 parsed clinical table",
    subtype_rule = paste(unique(na.omit(primary$subtype_mapping_rule)), collapse = ";"),
    expression_type = expression_note,
    annotation = annotation_note,
    min_signature_coverage = coverage_min,
    expected_sample_count_pass = expected_ok,
    endpoint_mapping_pass = endpoint_ok,
    subtype_mapping_pass = subtype_ok,
    interpretation_lock = ifelse(
      any(models$model_name == "axis_plus_proliferation" &
            models$variable == "scale(APC_CXCL9_axis)" &
            models$p_value >= 0.05, na.rm = TRUE),
      "supportive but not proliferation-independent evidence",
      "supportive and proliferation-independent pending manuscript review"
    ),
    notes = ifelse(expected_ok && endpoint_ok && subtype_ok,
                   "Mapping audit passed; use cautious supportive/proliferation-sensitive interpretation.",
                   "Mapping ambiguity detected; demote to exploratory supplement.")
  )
  readr::write_tsv(audit, file.path(results_dir, "GSE163882_mapping_audit.tsv"))

  p1 <- ggplot(primary, aes(x = factor(pCR_binary, levels = c(0, 1), labels = c("RD", "pCR")),
                            y = APC_CXCL9_axis, fill = factor(pCR_binary))) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.6) +
    geom_point(position = position_jitter(width = 0.12), size = 1.6, alpha = 0.75) +
    scale_fill_manual(values = c("0" = "#9CA3AF", "1" = "#2B6CB0")) +
    labs(x = NULL, y = "APC/CXCL9/HLA-II axis score") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", panel.grid.minor = element_blank())
  ggsave(file.path(figures_dir, "GSE163882_APC_axis_by_pCR.pdf"), p1, width = 4.8, height = 4.0, device = "pdf")

  p2 <- ggplot(primary, aes(x = Proliferation, y = APC_CXCL9_axis, color = factor(pCR_binary))) +
    geom_point(alpha = 0.78, size = 1.8) +
    geom_smooth(method = "lm", se = FALSE, color = "grey35", linewidth = 0.5) +
    scale_color_manual(values = c("0" = "#9CA3AF", "1" = "#2B6CB0"), labels = c("RD", "pCR")) +
    labs(x = "Proliferation score", y = "APC/CXCL9/HLA-II axis score", color = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave(file.path(figures_dir, "GSE163882_axis_vs_proliferation_scatter.pdf"), p2, width = 5.2, height = 4.2, device = "pdf")

  axis_row <- models %>% filter(model_name == "axis_only", variable == "scale(APC_CXCL9_axis)") %>% slice(1)
  adj_row <- models %>% filter(model_name == "axis_plus_proliferation", variable == "scale(APC_CXCL9_axis)") %>% slice(1)
  prolif_row <- models %>% filter(model_name == "axis_plus_proliferation", variable == "scale(Proliferation)") %>% slice(1)
  lines <- c(
    "# GSE163882 Mapping Audit Report",
    "",
    paste("Generated:", Sys.time()),
    "",
    "GSE163882 is the only statistically positive non-I-SPY2 external pCR cohort in Round 4, so Round 5 audits its clinical and expression mapping before manuscript use.",
    "",
    "## Mapping QC",
    "",
    paste(capture.output(print(audit, width = Inf)), collapse = "\n"),
    "",
    "## Clinical Variable Distribution",
    "",
    paste(capture.output(print(clinical_dist, n = Inf)), collapse = "\n"),
    "",
    "## Final Axis Models",
    "",
    paste0("- Axis-only OR=", fmt_num(axis_row$estimate), " (95% CI ", fmt_num(axis_row$conf_low), "-", fmt_num(axis_row$conf_high), "), p=", fmt_num(axis_row$p_value), "."),
    paste0("- Proliferation-adjusted axis OR=", fmt_num(adj_row$estimate), " (95% CI ", fmt_num(adj_row$conf_low), "-", fmt_num(adj_row$conf_high), "), p=", fmt_num(adj_row$p_value), "."),
    paste0("- Proliferation term OR=", fmt_num(prolif_row$estimate), ", p=", fmt_num(prolif_row$p_value), "."),
    "",
    "## Locked Interpretation",
    "",
    paste0("GSE163882 interpretation: `", audit$interpretation_lock, "`."),
    "Use this cohort as supportive evidence only; do not let it override GSE25066 as the technically validated negative chemotherapy-context cohort."
  )
  writeLines(lines, file.path(reports_dir, "GSE163882_mapping_audit_report.md"))
  message("GSE163882 mapping audit complete.")
}

main()
