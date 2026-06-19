#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
  library(yaml)
})

root <- getwd()
cfg <- yaml::read_yaml(file.path(root, "config", "project_config.yaml"))

as_binary_01 <- function(x) {
  x_chr <- trimws(as.character(x))
  dplyr::case_when(
    x_chr %in% c("1", "1.0", "TRUE", "True", "true", "yes", "Yes", "positive", "Positive") ~ 1,
    x_chr %in% c("0", "0.0", "FALSE", "False", "false", "no", "No", "negative", "Negative") ~ 0,
    TRUE ~ NA_real_
  )
}

find_required_col <- function(dat, candidates, label) {
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit) == 0) {
    stop("Could not find required ", label, " column. Candidates: ", paste(candidates, collapse = ", "))
  }
  hit[1]
}

normalize_treatment_arm <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x <- stringr::str_replace_all(x, "AMG-386", "AMG 386")
  x
}

map_ispy2_clinical <- function(pheno) {
  hr_col <- find_required_col(pheno, c("HR", "hr", "Hormone_receptor_status", "hormone_receptor_status"), "HR")
  her2_col <- find_required_col(pheno, c("HER2", "her2", "HER2_status", "her2_status"), "HER2")
  pcr_col <- find_required_col(pheno, c("pCR", "pcr", "PCR", "pathological_complete_response"), "pCR")
  arm_col <- find_required_col(pheno, c("Arm", "arm", "Treatment_arm", "treatment_arm"), "Arm")
  mp_col <- intersect(c("MP", "mp", "MammaPrint", "mammaprint"), names(pheno))
  mp_col <- if (length(mp_col) > 0) mp_col[1] else NA_character_

  mapped <- pheno %>%
    transmute(
      sample_id = sample_id,
      geo_accession = geo_accession,
      sample_title = sample_title,
      HR_raw = .data[[hr_col]],
      HER2_raw = .data[[her2_col]],
      pCR_raw = .data[[pcr_col]],
      MP_raw = if (!is.na(mp_col)) .data[[mp_col]] else NA_character_,
      treatment_arm_raw = as.character(.data[[arm_col]]),
      treatment_arm = normalize_treatment_arm(.data[[arm_col]])
    ) %>%
    mutate(
      HR_binary = as_binary_01(HR_raw),
      HER2_binary = as_binary_01(HER2_raw),
      pCR_binary = as_binary_01(pCR_raw),
      MP_binary = as_binary_01(MP_raw),
      HR_negative = HR_binary == 0,
      HER2_negative = HER2_binary == 0,
      tnbc_binary = HR_negative & HER2_negative,
      tnbc_definition = "GSE194040 primary proxy: HR == 0 and HER2 == 0",
      pCR_definition = "GSE194040 GEO annotation: pCR == 1 complete response, pCR == 0 failed complete response",
      notes = dplyr::case_when(
        is.na(HR_binary) | is.na(HER2_binary) ~ "TNBC proxy unavailable because HR or HER2 is missing",
        is.na(pCR_binary) ~ "pCR unavailable",
        TRUE ~ ""
      )
    )
  mapped
}

calc_auc <- function(y, pred) {
  ok <- !is.na(y) & !is.na(pred)
  y <- y[ok]
  pred <- pred[ok]
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  ranks <- rank(pred, ties.method = "average")
  (sum(ranks[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

roc_points <- function(y, pred) {
  ok <- !is.na(y) & !is.na(pred)
  y <- y[ok]
  pred <- pred[ok]
  thresholds <- sort(unique(pred), decreasing = TRUE)
  pts <- lapply(thresholds, function(th) {
    call <- pred >= th
    tp <- sum(call & y == 1)
    fp <- sum(call & y == 0)
    fn <- sum(!call & y == 1)
    tn <- sum(!call & y == 0)
    tibble(
      threshold = th,
      sensitivity = ifelse(tp + fn == 0, NA_real_, tp / (tp + fn)),
      specificity = ifelse(tn + fp == 0, NA_real_, tn / (tn + fp))
    )
  })
  bind_rows(tibble(threshold = Inf, sensitivity = 0, specificity = 1), bind_rows(pts), tibble(threshold = -Inf, sensitivity = 1, specificity = 0))
}

tidy_or_table <- function(model, model_name, dataset, endpoint, n_total, n_event, covariates) {
  broom::tidy(model) %>%
    mutate(
      dataset = dataset,
      analysis_subset = "HR0_HER20_TNBC_proxy",
      endpoint = endpoint,
      n_total = n_total,
      n_event_or_pCR = n_event,
      model_name = model_name,
      variable = term,
      effect_type = "OR",
      estimate = exp(estimate),
      conf_low = exp(log(estimate) - 1.96 * std.error),
      conf_high = exp(log(estimate) + 1.96 * std.error),
      p_value = p.value,
      covariates = covariates,
      notes = "Wald CI; treatment-arm coefficients retained for transparency.",
      .keep = "unused"
    ) %>%
    select(dataset, analysis_subset, endpoint, n_total, n_event_or_pCR, model_name, variable,
           effect_type, estimate, conf_low, conf_high, p_value, covariates, notes)
}

fit_logistic_models <- function(dat) {
  dat <- dat %>%
    filter(!is.na(pCR_binary), !is.na(TNBC_BIR_favorable_score)) %>%
    mutate(treatment_arm = factor(treatment_arm))
  if (nrow(dat) < 20) stop("Too few samples for pCR model after filtering: ", nrow(dat))
  if (sum(dat$pCR_binary == 1) < 5 || sum(dat$pCR_binary == 0) < 5) {
    stop("Too few pCR or residual-disease events for logistic modeling.")
  }

  arm_levels <- nlevels(dat$treatment_arm)
  use_arm <- arm_levels > 1 && nrow(dat) / arm_levels >= 5
  if (!use_arm) {
    warning("Treatment arm is not usable; fitting intercept-only baseline.")
  }

  f0 <- if (use_arm) pCR_binary ~ treatment_arm else pCR_binary ~ 1
  f1 <- if (use_arm) pCR_binary ~ treatment_arm + scale(TNBC_BIR_favorable_score) else pCR_binary ~ scale(TNBC_BIR_favorable_score)
  f2 <- if (use_arm) pCR_binary ~ treatment_arm + scale(Basoluminal_score) + scale(Immune_remodeling_score) + scale(Suppressive_myeloid_score) else pCR_binary ~ scale(Basoluminal_score) + scale(Immune_remodeling_score) + scale(Suppressive_myeloid_score)

  mods <- list(
    Model0_treatment_arm = glm(f0, family = binomial(), data = dat),
    Model1_treatment_arm_BIR = glm(f1, family = binomial(), data = dat),
    Model2_treatment_arm_components = glm(f2, family = binomial(), data = dat)
  )

  covariates <- c(
    Model0_treatment_arm = if (use_arm) "treatment_arm" else "none",
    Model1_treatment_arm_BIR = if (use_arm) "treatment_arm + TNBC_BIR_favorable_score" else "TNBC_BIR_favorable_score",
    Model2_treatment_arm_components = if (use_arm) "treatment_arm + Basoluminal_score + Immune_remodeling_score + Suppressive_myeloid_score" else "Basoluminal_score + Immune_remodeling_score + Suppressive_myeloid_score"
  )

  model_tbl <- purrr::imap_dfr(mods, function(m, nm) {
    tidy_or_table(
      m,
      model_name = nm,
      dataset = "GSE194040",
      endpoint = "pCR",
      n_total = nobs(m),
      n_event = sum(dat$pCR_binary == 1),
      covariates = covariates[[nm]]
    )
  })

  pred_tbl <- purrr::imap_dfr(mods, function(m, nm) {
    tibble(sample_id = dat$sample_id, model_name = nm, pCR_binary = dat$pCR_binary, pred = as.numeric(predict(m, type = "response")))
  })

  auc_tbl <- pred_tbl %>%
    group_by(model_name) %>%
    summarise(
      auc = calc_auc(pCR_binary, pred),
      n = n(),
      n_pCR = sum(pCR_binary == 1),
      n_residual = sum(pCR_binary == 0),
      .groups = "drop"
    )

  lrt_tbl <- tibble(
    comparison = c("Model1_vs_Model0", "Model2_vs_Model0"),
    p_value_lrt = c(
      anova(mods$Model0_treatment_arm, mods$Model1_treatment_arm_BIR, test = "Chisq")$`Pr(>Chi)`[2],
      anova(mods$Model0_treatment_arm, mods$Model2_treatment_arm_components, test = "Chisq")$`Pr(>Chi)`[2]
    )
  )

  list(models = mods, model_tbl = model_tbl, pred_tbl = pred_tbl, auc_tbl = auc_tbl, lrt_tbl = lrt_tbl, model_data = dat)
}

make_figures <- function(dat_tnbc, pred_tbl, auc_tbl) {
  fig_dir <- file.path(root, cfg$paths$figures)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  p1 <- dat_tnbc %>%
    mutate(pCR_label = if_else(pCR_binary == 1, "pCR", "Residual disease")) %>%
    ggplot(aes(x = pCR_label, y = TNBC_BIR_favorable_score, fill = pCR_label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55) +
    geom_jitter(width = 0.15, alpha = 0.65, size = 1.6) +
    labs(x = NULL, y = "TNBC-BIR favorable score", title = "GSE194040 TNBC: BIR score by pCR") +
    theme_bw() +
    theme(legend.position = "none")
  ggsave(file.path(fig_dir, "GSE194040_scores_by_pCR.pdf"), p1, width = 5.2, height = 4.2)

  roc_tbl <- pred_tbl %>%
    group_by(model_name) %>%
    group_modify(~ roc_points(.x$pCR_binary, .x$pred)) %>%
    ungroup() %>%
    left_join(auc_tbl, by = "model_name") %>%
    mutate(model_label = paste0(model_name, " (AUC=", sprintf("%.3f", auc), ")"))

  p2 <- roc_tbl %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, color = model_label)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray55") +
    geom_step(linewidth = 0.9) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(x = "1 - specificity", y = "Sensitivity", color = NULL, title = "GSE194040 TNBC pCR ROC") +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "GSE194040_ROC_models.pdf"), p2, width = 6.2, height = 5.0)

  quadrant_tbl <- dat_tnbc %>%
    group_by(BIR_quadrant) %>%
    summarise(n = n(), n_pCR = sum(pCR_binary == 1), pCR_rate = mean(pCR_binary == 1), .groups = "drop")

  p3 <- quadrant_tbl %>%
    ggplot(aes(x = reorder(BIR_quadrant, pCR_rate), y = pCR_rate, fill = BIR_quadrant)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = paste0("n=", n, ", pCR=", n_pCR)), hjust = -0.04, size = 3.2) +
    coord_flip() +
    scale_y_continuous(labels = percent_format(), limits = c(0, max(0.05, min(1, max(quadrant_tbl$pCR_rate) + 0.16)))) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = NULL, y = "pCR rate", title = "GSE194040 TNBC pCR rate by TNBC-BIR quadrant")
  ggsave(file.path(fig_dir, "GSE194040_quadrant_pCR_rate.pdf"), p3, width = 7.2, height = 4.6)

  calibration_tbl <- pred_tbl %>%
    filter(model_name == "Model1_treatment_arm_BIR") %>%
    mutate(bin = ntile(pred, 5)) %>%
    group_by(bin) %>%
    summarise(mean_pred = mean(pred), observed_pCR = mean(pCR_binary == 1), n = n(), .groups = "drop")
  p4 <- calibration_tbl %>%
    ggplot(aes(x = mean_pred, y = observed_pCR)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray55") +
    geom_point(size = 2.5) +
    geom_line() +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_bw() +
    labs(x = "Mean predicted pCR probability", y = "Observed pCR rate", title = "Model1 calibration by quintile")
  ggsave(file.path(fig_dir, "GSE194040_calibration_Model1.pdf"), p4, width = 5.2, height = 4.4)

  list(roc_tbl = roc_tbl, quadrant_tbl = quadrant_tbl, calibration_tbl = calibration_tbl)
}

main <- function() {
  processed_dir <- file.path(root, cfg$paths$processed_data)
  results_dir <- file.path(root, cfg$paths$results)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  scored <- readRDS(file.path(processed_dir, "GSE194040_scored_bulk.rds"))
  scores <- scored$scores
  pheno <- scored$pheno
  clin <- map_ispy2_clinical(pheno)
  readr::write_tsv(clin, file.path(processed_dir, "GSE194040_clinical.tsv"))

  dat_all <- scores %>% left_join(clin, by = "sample_id")
  readr::write_tsv(dat_all, file.path(results_dir, "GSE194040_all_scores_clinical.tsv"))

  dat_tnbc <- dat_all %>% filter(tnbc_binary %in% TRUE)
  readr::write_tsv(dat_tnbc, file.path(results_dir, "GSE194040_TNBC_scores_clinical.tsv"))

  arm_qc <- dat_tnbc %>%
    group_by(treatment_arm) %>%
    summarise(
      n = n(),
      n_pCR = sum(pCR_binary == 1, na.rm = TRUE),
      n_residual = sum(pCR_binary == 0, na.rm = TRUE),
      raw_arm_labels = paste(sort(unique(treatment_arm_raw)), collapse = " | "),
      .groups = "drop"
    ) %>%
    arrange(treatment_arm)
  readr::write_tsv(arm_qc, file.path(results_dir, "GSE194040_treatment_arm_qc.tsv"))

  count_flow <- tibble(
    dataset = "GSE194040",
    step = c(
      "signature_scored_samples",
      "clinical_samples",
      "score_clinical_joined_samples",
      "HR_mapped_samples",
      "HER2_mapped_samples",
      "pCR_mapped_samples",
      "TNBC_proxy_HR0_HER20_samples",
      "TNBC_proxy_with_pCR_and_BIR_complete"
    ),
    n = c(
      nrow(scores),
      nrow(clin),
      nrow(dat_all),
      sum(!is.na(dat_all$HR_binary)),
      sum(!is.na(dat_all$HER2_binary)),
      sum(!is.na(dat_all$pCR_binary)),
      nrow(dat_tnbc),
      sum(dat_tnbc$tnbc_binary %in% TRUE & !is.na(dat_tnbc$pCR_binary) & !is.na(dat_tnbc$TNBC_BIR_favorable_score))
    ),
    notes = c(
      "Rows in signature score table.",
      "Rows in parsed and mapped GEO clinical table.",
      "Left join by normalized I-SPY2 ResID/sample_id.",
      "HR parsed from GEO !Sample_characteristics_ch1.",
      "HER2 parsed from GEO !Sample_characteristics_ch1.",
      "pCR parsed from GEO !Sample_characteristics_ch1.",
      "Primary TNBC proxy per plan: HR == 0 and HER2 == 0.",
      "Final samples entering complete-case Model0/Model1/Model2 pCR analyses."
    )
  )
  readr::write_tsv(count_flow, file.path(results_dir, "GSE194040_sample_count_flow.tsv"))

  if (nrow(dat_tnbc) < 20) {
    stop("GSE194040 TNBC proxy subset has fewer than 20 samples after HR/HER2 mapping.")
  }

  res <- fit_logistic_models(dat_tnbc)
  readr::write_tsv(res$model_tbl, file.path(results_dir, "GSE194040_pCR_models.tsv"))
  readr::write_tsv(res$model_tbl, file.path(results_dir, "GSE194040_pCR_model_coefficients.tsv"))
  readr::write_tsv(res$auc_tbl, file.path(results_dir, "GSE194040_pCR_model_auc.tsv"))
  readr::write_tsv(res$lrt_tbl, file.path(results_dir, "GSE194040_pCR_model_lrt.tsv"))
  readr::write_tsv(res$pred_tbl, file.path(results_dir, "GSE194040_pCR_model_predictions.tsv"))

  direction_tbl <- res$model_data %>%
    mutate(pCR_label = if_else(pCR_binary == 1, "pCR", "Residual disease")) %>%
    group_by(pCR_label) %>%
    summarise(
      n = n(),
      mean_BIR = mean(TNBC_BIR_favorable_score, na.rm = TRUE),
      mean_Immune = mean(Immune_remodeling_score, na.rm = TRUE),
      mean_Basoluminal = mean(Basoluminal_score, na.rm = TRUE),
      mean_Suppressive_myeloid = mean(Suppressive_myeloid_score, na.rm = TRUE),
      .groups = "drop"
    )
  readr::write_tsv(direction_tbl, file.path(results_dir, "GSE194040_score_direction_by_pCR.tsv"))

  fig_tbls <- make_figures(res$model_data, res$pred_tbl, res$auc_tbl)
  readr::write_tsv(fig_tbls$roc_tbl, file.path(results_dir, "GSE194040_ROC_curve_points.tsv"))
  readr::write_tsv(fig_tbls$quadrant_tbl, file.path(results_dir, "GSE194040_quadrant_pCR_rate.tsv"))
  readr::write_tsv(fig_tbls$calibration_tbl, file.path(results_dir, "GSE194040_calibration_Model1.tsv"))

  saveRDS(res, file.path(processed_dir, "GSE194040_pCR_model_results.rds"))
  message("GSE194040 pCR modeling complete.")
}

main()
