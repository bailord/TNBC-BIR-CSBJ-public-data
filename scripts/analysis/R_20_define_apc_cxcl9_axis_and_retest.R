#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(survival)
  library(tibble)
  library(tidyr)
  library(yaml)
})

set.seed(42)
root <- getwd()
cfg <- yaml::read_yaml(file.path(root, "config", "project_config.yaml"))
results_dir <- file.path(root, cfg$paths$results)
fig_dir <- file.path(root, cfg$paths$figures)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

axis_components <- c("HLAII", "CXCL9_FOLR2_macrophage", "DC_HLA_antigen_presentation")

add_axis <- function(dat, group_vars = NULL) {
  calc <- function(x) {
    x %>%
      mutate(
        APC_CXCL9_axis = rowMeans(cbind(
          as.numeric(scale(HLAII)),
          as.numeric(scale(CXCL9_FOLR2_macrophage)),
          as.numeric(scale(DC_HLA_antigen_presentation))
        ), na.rm = TRUE)
      )
  }
  if (is.null(group_vars) || length(group_vars) == 0) return(calc(dat))
  dat %>% group_by(across(all_of(group_vars))) %>% group_modify(~ calc(.x)) %>% ungroup()
}

empty_model <- function(dataset, subset, endpoint, model_name, effect_type, notes) {
  tibble(dataset = dataset, analysis_subset = subset, endpoint = endpoint, n_total = NA_integer_, n_event_or_pCR = NA_integer_, model_name = model_name, variable = NA_character_, effect_type = effect_type, estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_, covariates = NA_character_, notes = notes)
}

fit_logistic <- function(dat, dataset, subset, model_name, formula, cohort_factor = NULL) {
  vars <- all.vars(formula)
  md <- dat %>% select(any_of(unique(c("sample_id", vars)))) %>% filter(if_all(all_of(vars), ~ !is.na(.)))
  if (!is.null(cohort_factor) && cohort_factor %in% names(md)) md[[cohort_factor]] <- factor(md[[cohort_factor]])
  n_event <- sum(md$pCR_binary == 1)
  n_rd <- sum(md$pCR_binary == 0)
  if (nrow(md) < 20 || n_event < 5 || n_rd < 5) return(empty_model(dataset, subset, "pCR", model_name, "OR", "Insufficient complete cases/events."))
  fit <- glm(formula, family = binomial(), data = md)
  broom::tidy(fit) %>%
    mutate(beta = estimate) %>%
    transmute(dataset = dataset, analysis_subset = subset, endpoint = "pCR", n_total = nrow(md), n_event_or_pCR = n_event, model_name = model_name, variable = term, effect_type = "OR", estimate = exp(beta), conf_low = exp(beta - 1.96 * std.error), conf_high = exp(beta + 1.96 * std.error), p_value = p.value, covariates = paste(deparse(formula), collapse = " "), notes = "")
}

fit_cox <- function(dat, dataset, subset, endpoint, model_name, formula, event_col) {
  vars <- all.vars(formula)
  md <- dat %>% select(any_of(unique(c("sample_id", vars)))) %>% filter(if_all(all_of(vars), ~ !is.na(.)))
  if ("cohort" %in% names(md)) md$cohort <- factor(md$cohort)
  n_event <- sum(md[[event_col]] == 1)
  if (nrow(md) < 30 || n_event < 10) return(empty_model(dataset, subset, endpoint, model_name, "HR", "Insufficient complete cases/events."))
  fit <- survival::coxph(formula, data = md)
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    transmute(dataset = dataset, analysis_subset = subset, endpoint = endpoint, n_total = nrow(md), n_event_or_pCR = n_event, model_name = model_name, variable = term, effect_type = "HR", estimate = estimate, conf_low = conf.low, conf_high = conf.high, p_value = p.value, covariates = paste(deparse(formula), collapse = " "), notes = "")
}

axis_correlations <- function(dat, dataset, subset) {
  dat2 <- dat %>% add_axis()
  vars <- intersect(c("TNBC_BIR_favorable_score", "Immune_remodeling_score", "HLAII", "Proliferation", "Cytolytic", "Basoluminal_score", "Suppressive_myeloid_score"), names(dat2))
  bind_rows(lapply(vars, function(v) {
    tibble(dataset = dataset, analysis_subset = subset, score = v, n = sum(complete.cases(dat2[, c("APC_CXCL9_axis", v)])), pearson = cor(dat2$APC_CXCL9_axis, dat2[[v]], use = "complete.obs", method = "pearson"), spearman = cor(dat2$APC_CXCL9_axis, dat2[[v]], use = "complete.obs", method = "spearman"))
  }))
}

main <- function() {
  g1940 <- readr::read_tsv(file.path(results_dir, "GSE194040_TNBC_scores_clinical.tsv"), show_col_types = FALSE) %>% add_axis()
  g250 <- readr::read_tsv(file.path(results_dir, "GSE25066_TNBC_scores_clinical.tsv"), show_col_types = FALSE) %>% add_axis(group_vars = "cohort")
  g250_rma <- if (file.exists(file.path(results_dir, "GSE25066_RMA_scores_clinical.tsv"))) {
    readr::read_tsv(file.path(results_dir, "GSE25066_RMA_scores_clinical.tsv"), show_col_types = FALSE) %>% add_axis(group_vars = "cohort")
  } else tibble()
  metabric <- readr::read_tsv(file.path(results_dir, "METABRIC_scores_clinical.tsv"), show_col_types = FALSE) %>% filter(clinical_tnbc %in% TRUE) %>% add_axis()

  g1940_models <- bind_rows(
    fit_logistic(g1940, "GSE194040", "primary_TNBC", "axis_treatment_adjusted", pCR_binary ~ treatment_arm + scale(APC_CXCL9_axis), "treatment_arm"),
    fit_logistic(g1940, "GSE194040", "primary_TNBC", "axis_plus_prolif_cytolytic", pCR_binary ~ treatment_arm + scale(APC_CXCL9_axis) + scale(Proliferation) + scale(Cytolytic), "treatment_arm"),
    bind_rows(lapply(sort(unique(g1940$treatment_arm)), function(arm) {
      sub <- g1940 %>% filter(treatment_arm == arm)
      fit_logistic(sub, "GSE194040", paste0("arm_", arm), "arm_specific_axis", pCR_binary ~ scale(APC_CXCL9_axis))
    }))
  )
  readr::write_tsv(g1940_models, file.path(results_dir, "GSE194040_APC_CXCL9_axis_models.tsv"))

  g250_primary <- g250 %>% filter(primary_tnbc %in% TRUE)
  g250_models <- bind_rows(
    fit_logistic(g250_primary, "GSE25066_processed", "primary_TNBC", "axis_cohort_adjusted", pCR_binary ~ cohort + scale(APC_CXCL9_axis), "cohort"),
    fit_logistic(g250_primary, "GSE25066_processed", "primary_TNBC", "axis_plus_proliferation", pCR_binary ~ cohort + scale(APC_CXCL9_axis) + scale(Proliferation), "cohort"),
    fit_cox(g250_primary, "GSE25066_processed", "primary_TNBC", "DRFS", "axis_cohort_adjusted", Surv(drfs_time_years, drfs_event) ~ cohort + scale(APC_CXCL9_axis), "drfs_event")
  )
  if (nrow(g250_rma) > 0) {
    rma_primary <- g250_rma %>% filter(primary_tnbc %in% TRUE)
    g250_models <- bind_rows(g250_models,
      fit_logistic(rma_primary, "GSE25066_RMA", "primary_TNBC", "axis_cohort_adjusted", pCR_binary ~ cohort + scale(APC_CXCL9_axis), "cohort"),
      fit_logistic(rma_primary, "GSE25066_RMA", "primary_TNBC", "axis_plus_proliferation", pCR_binary ~ cohort + scale(APC_CXCL9_axis) + scale(Proliferation), "cohort"),
      fit_cox(rma_primary, "GSE25066_RMA", "primary_TNBC", "DRFS", "axis_cohort_adjusted", Surv(drfs_time_years, drfs_event) ~ cohort + scale(APC_CXCL9_axis), "drfs_event")
    )
  }
  readr::write_tsv(g250_models, file.path(results_dir, "GSE25066_APC_CXCL9_axis_models.tsv"))

  metabric_models <- bind_rows(
    fit_cox(metabric, "METABRIC", "clinical_TNBC", "RFS", "axis_unadjusted", Surv(rfs_time_months, rfs_event) ~ scale(APC_CXCL9_axis), "rfs_event"),
    fit_cox(metabric, "METABRIC", "clinical_TNBC", "RFS", "axis_age_adjusted", Surv(rfs_time_months, rfs_event) ~ scale(APC_CXCL9_axis) + age, "rfs_event"),
    fit_cox(metabric, "METABRIC", "clinical_TNBC", "RFS", "axis_age_stage_complete_case", Surv(rfs_time_months, rfs_event) ~ scale(APC_CXCL9_axis) + age + tumor_stage, "rfs_event")
  )
  readr::write_tsv(metabric_models, file.path(results_dir, "METABRIC_APC_CXCL9_axis_survival.tsv"))

  cor_tbl <- bind_rows(
    axis_correlations(g1940, "GSE194040", "primary_TNBC"),
    axis_correlations(g250, "GSE25066_processed", "all_processed"),
    if (nrow(g250_rma) > 0) axis_correlations(g250_rma, "GSE25066_RMA", "all_RMA") else tibble(),
    axis_correlations(metabric, "METABRIC", "clinical_TNBC")
  )
  readr::write_tsv(cor_tbl, file.path(results_dir, "APC_CXCL9_axis_correlation_with_existing_scores.tsv"))

  summary_effects <- bind_rows(g1940_models, g250_models, metabric_models) %>%
    filter(variable == "scale(APC_CXCL9_axis)") %>%
    mutate(label = paste(dataset, endpoint, model_name, sep = " | "))
  if (nrow(summary_effects) > 0) {
    p <- ggplot(summary_effects, aes(x = estimate, y = reorder(label, estimate), color = endpoint)) +
      geom_vline(xintercept = 1, linetype = 2, color = "gray55") +
      geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0.15, orientation = "y") +
      geom_point(size = 2) +
      scale_x_log10() +
      theme_bw() +
      labs(x = "Effect per 1 SD APC/CXCL9/HLA-II axis", y = NULL, title = "APC/CXCL9/HLA-II axis cross-dataset summary")
    ggsave(file.path(fig_dir, "APC_CXCL9_axis_cross_dataset_summary.pdf"), p, width = 9, height = max(5, 0.3 * nrow(summary_effects)))
  }
  message("APC/CXCL9/HLA-II axis retest complete.")
}

main()
