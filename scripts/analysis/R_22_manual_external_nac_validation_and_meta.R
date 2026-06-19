#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(metafor)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(yaml)
})

set.seed(42)
root <- getwd()
cfg <- yaml::read_yaml(file.path(root, "config", "project_config.yaml"))
results_dir <- file.path(root, cfg$paths$results)
processed_dir <- file.path(root, cfg$paths$processed_data)
reports_dir <- file.path(root, "reports")
figures_dir <- file.path(root, "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

model_empty <- tibble(
  dataset = character(),
  analysis_subset = character(),
  endpoint = character(),
  n_total = integer(),
  n_event_or_pCR = integer(),
  n_non_event_or_RD = integer(),
  model_name = character(),
  variable = character(),
  effect_type = character(),
  estimate = double(),
  conf_low = double(),
  conf_high = double(),
  p_value = double(),
  log_effect = double(),
  se_log_effect = double(),
  covariates = character(),
  source_context = character(),
  notes = character()
)

meta_empty <- tibble(
  meta_set = character(),
  include_GSE194040 = logical(),
  include_GSE25066 = logical(),
  endpoint = character(),
  analysis_subset = character(),
  variable = character(),
  effect_type = character(),
  k = integer(),
  fixed_estimate = double(),
  fixed_conf_low = double(),
  fixed_conf_high = double(),
  fixed_p_value = double(),
  random_estimate = double(),
  random_conf_low = double(),
  random_conf_high = double(),
  random_p_value = double(),
  I2 = double(),
  Q_p_value = double(),
  included_datasets = character(),
  notes = character()
)

safe_scale <- function(x) {
  out <- as.numeric(scale(x))
  out[is.na(out)] <- NA_real_
  out
}

wald_from_fit <- function(fit, data, dataset, analysis_subset, model_name, term_map, covariates, notes) {
  sm <- coef(summary(fit))
  rows <- lapply(names(term_map), function(term) {
    if (!term %in% rownames(sm)) return(NULL)
    beta <- unname(coef(fit)[term])
    se <- unname(sm[term, "Std. Error"])
    tibble(
      dataset = dataset,
      analysis_subset = analysis_subset,
      endpoint = "pCR",
      n_total = nrow(data),
      n_event_or_pCR = sum(data$pCR_binary == 1, na.rm = TRUE),
      n_non_event_or_RD = sum(data$pCR_binary == 0, na.rm = TRUE),
      model_name = model_name,
      variable = term_map[[term]],
      effect_type = "OR_per_1SD",
      estimate = exp(beta),
      conf_low = exp(beta - 1.96 * se),
      conf_high = exp(beta + 1.96 * se),
      p_value = unname(sm[term, "Pr(>|z|)"]),
      log_effect = beta,
      se_log_effect = se,
      covariates = covariates,
      source_context = "manual_external_NAC",
      notes = notes
    )
  })
  bind_rows(rows)
}

fit_logistic_model <- function(df, dataset, analysis_subset, model_name, variables) {
  needed <- c("pCR_binary", variables)
  fit_df <- df %>%
    select(any_of(c(needed, "platform"))) %>%
    filter(if_all(all_of(needed), ~ !is.na(.x)))
  platform_adjust <- "platform" %in% names(fit_df) &&
    all(!is.na(fit_df$platform)) &&
    n_distinct(fit_df$platform) > 1
  if (platform_adjust) {
    fit_df <- fit_df %>% mutate(platform = factor(platform))
  } else {
    fit_df <- fit_df %>% select(-any_of("platform"))
  }
  n_pcr <- sum(fit_df$pCR_binary == 1, na.rm = TRUE)
  n_rd <- sum(fit_df$pCR_binary == 0, na.rm = TRUE)
  notes <- c()
  if (nrow(fit_df) < 25 || n_pcr < 8 || n_rd < 8) {
    return(model_empty %>% add_row(
      dataset = dataset, analysis_subset = analysis_subset, endpoint = "pCR",
      n_total = nrow(fit_df), n_event_or_pCR = n_pcr, n_non_event_or_RD = n_rd,
      model_name = model_name, variable = paste(variables, collapse = "+"),
      effect_type = "OR_per_1SD", estimate = NA_real_, conf_low = NA_real_,
      conf_high = NA_real_, p_value = NA_real_, log_effect = NA_real_,
      se_log_effect = NA_real_, covariates = "not_fit", source_context = "manual_external_NAC",
      notes = "Not fit: eligibility requires n>=25 and pCR/RD each >=8."
    ))
  }
  if (n_pcr < 15) notes <- c(notes, "underpowered: pCR events <15")
  if (platform_adjust) notes <- c(notes, "platform fixed effect included")
  for (v in variables) fit_df[[paste0("z_", v)]] <- safe_scale(fit_df[[v]])
  zvars <- paste0("z_", variables)
  formula_terms <- c(zvars, if (platform_adjust) "platform")
  formula <- as.formula(paste("pCR_binary ~", paste(formula_terms, collapse = " + ")))
  fit <- tryCatch(suppressWarnings(glm(formula, data = fit_df, family = binomial())), error = function(e) e)
  if (inherits(fit, "error")) {
    return(model_empty %>% add_row(
      dataset = dataset, analysis_subset = analysis_subset, endpoint = "pCR",
      n_total = nrow(fit_df), n_event_or_pCR = n_pcr, n_non_event_or_RD = n_rd,
      model_name = model_name, variable = paste(variables, collapse = "+"),
      effect_type = "OR_per_1SD", estimate = NA_real_, conf_low = NA_real_,
      conf_high = NA_real_, p_value = NA_real_, log_effect = NA_real_,
      se_log_effect = NA_real_, covariates = deparse(formula), source_context = "manual_external_NAC",
      notes = paste("glm failed:", conditionMessage(fit))
    ))
  }
  term_map <- stats::setNames(variables, zvars)
  wald_from_fit(
    fit, fit_df, dataset, analysis_subset, model_name, term_map,
    covariates = deparse(formula),
    notes = paste(notes, collapse = "; ")
  )
}

subset_specs <- list(
  primary_TNBC = function(x) x %>% filter(primary_tnbc, pCR_mapped),
  ER_HER2_negative_sensitivity = function(x) x %>% filter(er_her2_negative, pCR_mapped)
)

run_models <- function(scores, eligible_datasets) {
  rows <- list()
  for (dataset in eligible_datasets) {
    df0 <- scores %>% filter(dataset == !!dataset)
    for (subset_name in names(subset_specs)) {
      df <- subset_specs[[subset_name]](df0)
      rows[[length(rows) + 1]] <- fit_logistic_model(df, dataset, subset_name, "axis_only", "APC_CXCL9_axis")
      rows[[length(rows) + 1]] <- fit_logistic_model(df, dataset, subset_name, "axis_plus_proliferation", c("APC_CXCL9_axis", "Proliferation"))
      rows[[length(rows) + 1]] <- fit_logistic_model(df, dataset, subset_name, "immune_only", "Immune_remodeling_score")
      rows[[length(rows) + 1]] <- fit_logistic_model(df, dataset, subset_name, "proliferation_only", "Proliferation")
      rows[[length(rows) + 1]] <- fit_logistic_model(df, dataset, subset_name, "aggregate_BIR_descriptive", "TNBC_BIR_favorable_score")
    }
  }
  bind_rows(rows)
}

effect_from_axis_file <- function(path, dataset_pattern = NULL, model_name = NULL) {
  if (!file.exists(path)) return(tibble())
  tbl <- readr::read_tsv(path, show_col_types = FALSE)
  out <- tbl %>%
    filter(endpoint == "pCR", grepl("primary", analysis_subset, ignore.case = TRUE),
           variable == "scale(APC_CXCL9_axis)") %>%
    mutate(
      log_effect = log(estimate),
      se_log_effect = (log(conf_high) - log(conf_low)) / (2 * 1.96),
      n_non_event_or_RD = n_total - n_event_or_pCR,
      effect_type = "OR_per_1SD"
    )
  if (!is.null(dataset_pattern)) out <- out %>% filter(grepl(dataset_pattern, dataset))
  if (!is.null(model_name)) out <- out %>% filter(.data$model_name == !!model_name)
  out %>%
    transmute(
      dataset, analysis_subset, endpoint, n_total, n_event_or_pCR, n_non_event_or_RD,
      model_name, variable = "APC_CXCL9_axis", effect_type, estimate, conf_low, conf_high,
      p_value, log_effect, se_log_effect, covariates, source_context = "context_existing_external_validation",
      notes = ifelse(is.na(notes), "", notes)
    )
}

fit_meta <- function(effects, meta_set, include_gse194040, include_gse25066, notes) {
  effects <- effects %>%
    filter(is.finite(log_effect), is.finite(se_log_effect), se_log_effect > 0)
  if (nrow(effects) < 2) {
    return(meta_empty %>% add_row(
      meta_set = meta_set, include_GSE194040 = include_gse194040, include_GSE25066 = include_gse25066,
      endpoint = "pCR", analysis_subset = "primary_TNBC", variable = "APC_CXCL9_axis",
      effect_type = "OR_per_1SD", k = nrow(effects),
      fixed_estimate = NA_real_, fixed_conf_low = NA_real_, fixed_conf_high = NA_real_, fixed_p_value = NA_real_,
      random_estimate = NA_real_, random_conf_low = NA_real_, random_conf_high = NA_real_, random_p_value = NA_real_,
      I2 = NA_real_, Q_p_value = NA_real_, included_datasets = paste(effects$dataset, collapse = ";"),
      notes = paste(notes, "Meta-analysis not fit because k<2.")
    ))
  }
  fe <- metafor::rma.uni(yi = effects$log_effect, sei = effects$se_log_effect, method = "FE")
  re <- metafor::rma.uni(yi = effects$log_effect, sei = effects$se_log_effect, method = "REML")
  meta_empty %>% add_row(
    meta_set = meta_set,
    include_GSE194040 = include_gse194040,
    include_GSE25066 = include_gse25066,
    endpoint = "pCR",
    analysis_subset = "primary_TNBC",
    variable = "APC_CXCL9_axis",
    effect_type = "OR_per_1SD",
    k = nrow(effects),
    fixed_estimate = exp(as.numeric(fe$b)),
    fixed_conf_low = exp(fe$ci.lb),
    fixed_conf_high = exp(fe$ci.ub),
    fixed_p_value = fe$pval,
    random_estimate = exp(as.numeric(re$b)),
    random_conf_low = exp(re$ci.lb),
    random_conf_high = exp(re$ci.ub),
    random_p_value = re$pval,
    I2 = re$I2,
    Q_p_value = re$QEp,
    included_datasets = paste(effects$dataset, collapse = ";"),
    notes = notes
  )
}

make_forest_pdf <- function(effects, path) {
  if (!nrow(effects)) {
    pdf(path, width = 7, height = 4)
    plot.new()
    text(0.5, 0.5, "No eligible external NAC effects to plot.")
    dev.off()
    return(invisible(FALSE))
  }
  plot_df <- effects %>%
    mutate(
      dataset_label = paste0(dataset, " (n=", n_total, ", pCR=", n_event_or_pCR, ")"),
      context = case_when(
        dataset == "GSE194040" ~ "I-SPY2 context",
        grepl("GSE25066", dataset) ~ "GSE25066 negative context",
        TRUE ~ "manual external NAC"
      ),
      dataset_label = factor(dataset_label, levels = rev(dataset_label))
    )
  p <- ggplot(plot_df, aes(x = estimate, y = dataset_label, xmin = conf_low, xmax = conf_high, color = context)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4, color = "grey45") +
    geom_pointrange(linewidth = 0.55) +
    scale_x_log10() +
    labs(x = "OR per 1 SD APC/CXCL9/HLA-II axis", y = NULL, color = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave(path, p, width = 7.2, height = max(4.2, 0.38 * nrow(plot_df) + 1.2), device = "pdf")
  invisible(TRUE)
}

make_meta_pdf <- function(meta_tbl, path) {
  plot_df <- meta_tbl %>%
    filter(is.finite(random_estimate)) %>%
    transmute(
      meta_set,
      estimate = random_estimate,
      conf_low = random_conf_low,
      conf_high = random_conf_high,
      label = paste0(meta_set, " (k=", k, ")")
    )
  if (!nrow(plot_df)) {
    pdf(path, width = 7, height = 4)
    plot.new()
    text(0.5, 0.5, "No meta-analysis result to plot.")
    dev.off()
    return(invisible(FALSE))
  }
  plot_df$label <- factor(plot_df$label, levels = rev(plot_df$label))
  p <- ggplot(plot_df, aes(x = estimate, y = label, xmin = conf_low, xmax = conf_high)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4, color = "grey45") +
    geom_pointrange(color = "#2B6CB0", linewidth = 0.6) +
    scale_x_log10() +
    labs(x = "Random-effects pooled OR per 1 SD APC/CXCL9/HLA-II axis", y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())
  ggsave(path, p, width = 7.4, height = max(3.6, 0.45 * nrow(plot_df) + 1.1), device = "pdf")
  invisible(TRUE)
}

write_notes <- function(models, meta_tbl, qc) {
  primary_axis <- models %>%
    filter(analysis_subset == "primary_TNBC", model_name == "axis_only", variable == "APC_CXCL9_axis")
  support <- primary_axis %>%
    mutate(supports_direction = is.finite(estimate) & estimate > 1) %>%
    transmute(dataset, n_total, n_event_or_pCR, estimate, conf_low, conf_high, p_value, supports_direction, notes)
  lines <- c(
    "# Round 4 External NAC Manual Validation Notes",
    "",
    paste("Generated:", Sys.time()),
    "",
    "Manual external NAC modeling used only cohorts marked `eligible_for_external_validation_modeling` by `external_validation_external_dataset_mapping_qc.tsv`. No dataset was added or removed based on model direction.",
    "",
    "## Primary TNBC APC/CXCL9 Axis Results",
    "",
    paste(capture.output(print(support, n = Inf)), collapse = "\n"),
    "",
    "## Meta Summary",
    "",
    paste(capture.output(print(meta_tbl, n = Inf)), collapse = "\n"),
    "",
    "## Interpretation Guardrails",
    "",
    "- GSE25066 is shown as a decisive negative context cohort and is not treated as a rescue success.",
    "- GSE22226 and GSE32646 are retained when they pass the minimum rule, but rows with pCR events <15 are explicitly flagged underpowered.",
    "- Aggregate BIR appears only as descriptive/supplementary in this script."
  )
  writeLines(lines, file.path(reports_dir, "external_NAC_manual_validation_and_meta_report.md"))
}

main <- function() {
  mapping_path <- file.path(results_dir, "external_validation_external_dataset_mapping_qc.tsv")
  score_path <- file.path(processed_dir, "external_validation_external_NAC_scores_clinical.tsv")
  if (!file.exists(mapping_path)) stop("Run scripts/R/21_manual_external_nac_dataset_mapper.R first.")
  if (!file.exists(score_path)) stop("Missing processed external NAC score table.")
  qc <- readr::read_tsv(mapping_path, show_col_types = FALSE)
  scores <- readr::read_tsv(score_path, show_col_types = FALSE)
  eligible <- qc %>%
    filter(decision == "eligible_for_external_validation_modeling") %>%
    pull(dataset)

  if (!length(eligible)) {
    note <- "No manual external NAC cohort passed Round 4 eligibility after mapping; no external rescue model was run."
    readr::write_tsv(model_empty %>% add_row(
      dataset = NA_character_, analysis_subset = NA_character_, endpoint = "pCR",
      n_total = NA_integer_, n_event_or_pCR = NA_integer_, n_non_event_or_RD = NA_integer_,
      model_name = NA_character_, variable = NA_character_, effect_type = "OR_per_1SD",
      estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_,
      log_effect = NA_real_, se_log_effect = NA_real_, covariates = NA_character_,
      source_context = "manual_external_NAC", notes = note
    ), file.path(results_dir, "external_validation_external_NAC_models.tsv"))
    readr::write_tsv(meta_empty %>% add_row(
      meta_set = "manual_external_primary_TNBC", include_GSE194040 = FALSE, include_GSE25066 = FALSE,
      endpoint = "pCR", analysis_subset = "primary_TNBC", variable = "APC_CXCL9_axis",
      effect_type = "OR_per_1SD", k = 0L, notes = note
    ), file.path(results_dir, "external_validation_external_NAC_meta_summary.tsv"))
    make_forest_pdf(tibble(), file.path(figures_dir, "external_NAC_forest_APC_CXCL9_axis.pdf"))
    make_meta_pdf(tibble(), file.path(figures_dir, "external_NAC_pCR_meta_forest.pdf"))
    return(invisible(FALSE))
  }

  models <- run_models(scores, eligible)
  readr::write_tsv(models, file.path(results_dir, "external_validation_external_NAC_models.tsv"))

  manual_axis <- models %>%
    filter(
      analysis_subset == "primary_TNBC",
      model_name == "axis_only",
      variable == "APC_CXCL9_axis",
      is.finite(log_effect),
      is.finite(se_log_effect)
    )

  gse194040_axis <- effect_from_axis_file(
    file.path(results_dir, "GSE194040_APC_CXCL9_axis_models.tsv"),
    dataset_pattern = "^GSE194040$",
    model_name = "axis_treatment_adjusted"
  ) %>%
    filter(variable == "APC_CXCL9_axis")
  gse25066_axis <- effect_from_axis_file(
    file.path(results_dir, "GSE25066_APC_CXCL9_axis_models.tsv"),
    dataset_pattern = "^GSE25066_processed$",
    model_name = "axis_cohort_adjusted"
  ) %>%
    filter(variable == "APC_CXCL9_axis")

  meta_tbl <- bind_rows(
    fit_meta(manual_axis, "manual_external_primary_TNBC", FALSE, FALSE,
             "Non-I-SPY2 manual external NAC cohorts only."),
    fit_meta(bind_rows(manual_axis, gse194040_axis), "manual_plus_GSE194040_context", TRUE, FALSE,
             "Context-only pooled result including I-SPY2/GSE194040; not used as non-I-SPY2 validation."),
    fit_meta(bind_rows(manual_axis, gse25066_axis), "manual_plus_GSE25066_processed_context", FALSE, TRUE,
             "Context-only pooled result including decisive negative GSE25066 processed axis model.")
  )
  readr::write_tsv(meta_tbl, file.path(results_dir, "external_validation_external_NAC_meta_summary.tsv"))

  forest_effects <- bind_rows(
    manual_axis,
    gse194040_axis,
    gse25066_axis
  ) %>%
    mutate(
      dataset = recode(dataset, GSE25066_processed = "GSE25066_processed")
    )
  readr::write_tsv(forest_effects, file.path(results_dir, "external_validation_external_NAC_forest_effects.tsv"))
  make_forest_pdf(forest_effects, file.path(figures_dir, "external_NAC_forest_APC_CXCL9_axis.pdf"))
  make_meta_pdf(meta_tbl, file.path(figures_dir, "external_NAC_pCR_meta_forest.pdf"))
  write_notes(models, meta_tbl, qc)

  message("Round 4 external NAC manual validation and meta-analysis complete.")
  invisible(TRUE)
}

main()
