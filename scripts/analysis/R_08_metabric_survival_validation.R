#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(readr)
  library(stringr)
  library(survival)
  library(tibble)
  library(yaml)
})

source("scripts/R/utils_signature_scoring.R")

root <- getwd()
cfg <- yaml::read_yaml(file.path(root, "config", "project_config.yaml"))
public_root <- cfg$paths$public_data_root
metabric_dir <- file.path(public_root, "cBioPortal", "brca_metabric")
expr_path <- file.path(metabric_dir, "data_mrna_illumina_microarray.txt")
sample_path <- file.path(metabric_dir, "data_clinical_sample.txt")
patient_path <- file.path(metabric_dir, "data_clinical_patient.txt")
sig_tbl <- readr::read_tsv(file.path(root, "data/processed/TNBC_BIR_signatures_clean.tsv"), show_col_types = FALSE)

sha256_file <- function(path) {
  out <- tryCatch(system2("shasum", c("-a", "256", path), stdout = TRUE), error = function(e) NA_character_)
  if (length(out) == 0 || is.na(out[1])) return(NA_character_)
  strsplit(out[1], "\\s+")[[1]][1]
}

append_public_manifest <- function(paths) {
  manifest_path <- file.path(root, "results", "public_data_acquisition_manifest.tsv")
  rows <- tibble(
    dataset = "METABRIC",
    source_url = "Local cBioPortal public dataset already present under public_bioinf_databases/cBioPortal/brca_metabric",
    public_path = normalizePath(paths, mustWork = TRUE),
    sha256 = vapply(paths, sha256_file, character(1)),
    download_date = format(Sys.Date()),
    used_in_phase = "Phase5_METABRIC_survival_validation"
  )
  old <- if (file.exists(manifest_path)) {
    readr::read_tsv(manifest_path, show_col_types = FALSE) %>%
      mutate(download_date = as.character(download_date))
  } else {
    tibble()
  }
  bind_rows(old, rows) %>%
    distinct(dataset, public_path, .keep_all = TRUE) %>%
    arrange(dataset, basename(public_path)) %>%
    readr::write_tsv(manifest_path)
}

read_metabric_expression <- function(path) {
  message("Reading METABRIC expression: ", path)
  expr_tbl <- readr::read_tsv(path, show_col_types = FALSE, progress = TRUE)
  gene <- toupper(trimws(as.character(expr_tbl$Hugo_Symbol)))
  mat <- as.matrix(as.data.frame(lapply(expr_tbl[, -(1:2)], as.numeric), check.names = FALSE))
  rownames(mat) <- gene
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  mat
}

read_cbio_clinical <- function(path) {
  readr::read_tsv(path, comment = "#", show_col_types = FALSE)
}

event_from_status <- function(x) {
  dplyr::case_when(
    str_detect(as.character(x), "^1:") ~ 1,
    str_detect(as.character(x), "^0:") ~ 0,
    TRUE ~ NA_real_
  )
}

map_metabric_clinical <- function(sample_tbl, patient_tbl) {
  sample_tbl %>%
    left_join(patient_tbl, by = "PATIENT_ID") %>%
    transmute(
      sample_id = SAMPLE_ID,
      patient_id = PATIENT_ID,
      ER_status = ER_STATUS,
      PR_status = PR_STATUS,
      HER2_status = HER2_STATUS,
      ER_IHC = ER_IHC,
      HER2_SNP6 = HER2_SNP6,
      PAM50 = CLAUDIN_SUBTYPE,
      threegene = THREEGENE,
      age = suppressWarnings(as.numeric(AGE_AT_DIAGNOSIS)),
      tumor_stage = TUMOR_STAGE,
      rfs_time_months = suppressWarnings(as.numeric(RFS_MONTHS)),
      rfs_event = event_from_status(RFS_STATUS),
      os_time_months = suppressWarnings(as.numeric(OS_MONTHS)),
      os_event = event_from_status(OS_STATUS),
      ER_negative = str_detect(ER_status, regex("^Negative$", ignore_case = TRUE)),
      PR_negative = str_detect(PR_status, regex("^Negative$", ignore_case = TRUE)),
      HER2_negative = str_detect(HER2_status, regex("^Negative$", ignore_case = TRUE)),
      clinical_tnbc = ER_negative & PR_negative & HER2_negative,
      basal_like_proxy = str_detect(PAM50, regex("^Basal$|claudin-low", ignore_case = TRUE)) |
        str_detect(threegene, regex("ER-/HER2-", ignore_case = TRUE)),
      tnbc_definition = "Clinical TNBC: ER negative, PR negative, HER2 negative from cBioPortal sample table; basal-like proxy kept as sensitivity only."
    )
}

make_model_table <- function(fit, dat, dataset, subset_name, endpoint, model_name, covariates) {
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    transmute(
      dataset = dataset,
      analysis_subset = subset_name,
      endpoint = endpoint,
      n_total = nrow(dat),
      n_event_or_pCR = sum(dat$event == 1, na.rm = TRUE),
      model_name = model_name,
      variable = term,
      effect_type = "HR",
      estimate = estimate,
      conf_low = conf.low,
      conf_high = conf.high,
      p_value = p.value,
      covariates = covariates,
      notes = "Cox proportional hazards model; continuous scores scaled within METABRIC."
    )
}

fit_endpoint <- function(dat, subset_name, endpoint_name, time_col, event_col) {
  dat <- dat %>%
    transmute(
      sample_id, TNBC_BIR_favorable_score, Basoluminal_score,
      Immune_remodeling_score, Suppressive_myeloid_score,
      time = .data[[time_col]],
      event = .data[[event_col]]
    ) %>%
    filter(!is.na(time), !is.na(event), !is.na(TNBC_BIR_favorable_score))

  if (nrow(dat) < 30 || sum(dat$event == 1, na.rm = TRUE) < 10) {
    return(tibble(
      dataset = "METABRIC",
      analysis_subset = subset_name,
      endpoint = endpoint_name,
      n_total = nrow(dat),
      n_event_or_pCR = sum(dat$event == 1, na.rm = TRUE),
      model_name = NA_character_,
      variable = NA_character_,
      effect_type = "HR",
      estimate = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      covariates = NA_character_,
      notes = "Insufficient complete cases or events."
    ))
  }

  m1 <- coxph(Surv(time, event) ~ scale(TNBC_BIR_favorable_score), data = dat)
  m2 <- coxph(Surv(time, event) ~ scale(Basoluminal_score) + scale(Immune_remodeling_score) + scale(Suppressive_myeloid_score), data = dat)
  bind_rows(
    make_model_table(m1, dat, "METABRIC", subset_name, endpoint_name, "BIR_favorable_continuous", "TNBC_BIR_favorable_score"),
    make_model_table(m2, dat, "METABRIC", subset_name, endpoint_name, "components_continuous", "Basoluminal_score + Immune_remodeling_score + Suppressive_myeloid_score")
  )
}

plot_km <- function(dat, subset_name, endpoint_name, time_col, event_col, out_file) {
  dat <- dat %>%
    transmute(
      time = .data[[time_col]],
      event = .data[[event_col]],
      TNBC_BIR_favorable_score = TNBC_BIR_favorable_score
    ) %>%
    filter(!is.na(time), !is.na(event), !is.na(TNBC_BIR_favorable_score)) %>%
    mutate(
      BIR_tertile = ntile(TNBC_BIR_favorable_score, 3),
      BIR_group = factor(
        c("Low favorable", "Intermediate", "High favorable")[BIR_tertile],
        levels = c("Low favorable", "Intermediate", "High favorable")
      )
    )
  if (nrow(dat) < 30 || sum(dat$event == 1, na.rm = TRUE) < 10) return(invisible(FALSE))
  fit <- survfit(Surv(time, event) ~ BIR_group, data = dat)
  lr <- survdiff(Surv(time, event) ~ BIR_group, data = dat)
  p <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)

  pdf(out_file, width = 6.4, height = 5.2)
  plot(
    fit,
    col = c("#0072B2", "#666666", "#D55E00"),
    lwd = 2,
    xlab = "Months",
    ylab = paste(endpoint_name, "probability"),
    main = paste("METABRIC", subset_name, endpoint_name, "by BIR tertile")
  )
  legend("bottomleft", legend = levels(dat$BIR_group), col = c("#0072B2", "#666666", "#D55E00"), lwd = 2, bty = "n")
  text(x = max(dat$time, na.rm = TRUE) * 0.62, y = 0.12, labels = paste0("log-rank p=", signif(p, 3)))
  dev.off()
  invisible(TRUE)
}

main <- function() {
  stopifnot(file.exists(expr_path), file.exists(sample_path), file.exists(patient_path))
  dir.create(file.path(root, "data", "processed"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "figures"), recursive = TRUE, showWarnings = FALSE)
  append_public_manifest(c(expr_path, sample_path, patient_path))

  expr <- read_metabric_expression(expr_path)
  sample_tbl <- read_cbio_clinical(sample_path)
  patient_tbl <- read_cbio_clinical(patient_path)
  clin <- map_metabric_clinical(sample_tbl, patient_tbl)

  common <- intersect(colnames(expr), clin$sample_id)
  expr <- expr[, common, drop = FALSE]
  clin <- clin %>% filter(sample_id %in% common)

  res <- score_signatures(
    expr,
    sig_tbl = sig_tbl,
    min_coverage = cfg$score_settings$min_signature_coverage,
    duplicate_rule = cfg$score_settings$duplicate_gene_rule
  )
  scores <- assign_bir_quadrants(res$scores)
  coverage <- res$coverage %>% mutate(dataset = "METABRIC", .before = 1)
  dat <- scores %>% left_join(clin, by = c("sample_id" = "sample_id"))

  readr::write_tsv(clin, file.path(root, "data", "processed", "METABRIC_clinical.tsv"))
  readr::write_tsv(scores, file.path(root, "results", "METABRIC_TNBC_BIR_scores.tsv"))
  readr::write_tsv(coverage, file.path(root, "results", "METABRIC_signature_coverage.tsv"))
  readr::write_tsv(dat, file.path(root, "results", "METABRIC_scores_clinical.tsv"))
  saveRDS(list(expr = expr, clin = clin, scores = scores, coverage = coverage, scores_clinical = dat),
          file.path(root, "data", "processed", "METABRIC_scored_bulk.rds"))

  count_flow <- tibble(
    dataset = "METABRIC",
    step = c(
      "expression_samples",
      "clinical_samples",
      "expression_clinical_overlap",
      "clinical_TNBC_samples",
      "basal_like_proxy_samples",
      "clinical_TNBC_RFS_complete",
      "clinical_TNBC_OS_complete"
    ),
    n = c(
      ncol(expr),
      nrow(clin),
      nrow(dat),
      sum(dat$clinical_tnbc %in% TRUE, na.rm = TRUE),
      sum(dat$basal_like_proxy %in% TRUE, na.rm = TRUE),
      sum(dat$clinical_tnbc %in% TRUE & !is.na(dat$rfs_time_months) & !is.na(dat$rfs_event) & !is.na(dat$TNBC_BIR_favorable_score)),
      sum(dat$clinical_tnbc %in% TRUE & !is.na(dat$os_time_months) & !is.na(dat$os_event) & !is.na(dat$TNBC_BIR_favorable_score))
    ),
    notes = c(
      "Expression columns in cBioPortal METABRIC mRNA microarray matrix after duplicate-gene collapse during scoring.",
      "Clinical rows with matching expression sample IDs after sample/patient table join.",
      "Rows retained after expression-clinical sample intersection.",
      "Primary clinical TNBC: ER-, PR-, HER2- from sample-level cBioPortal annotations.",
      "Sensitivity proxy: PAM50 Basal or claudin-low and/or 3-gene ER-/HER2-.",
      "Complete clinical TNBC cases for RFS Cox model.",
      "Complete clinical TNBC cases for OS Cox model."
    )
  )
  readr::write_tsv(count_flow, file.path(root, "results", "METABRIC_sample_count_flow.tsv"))

  tnbc <- dat %>% filter(clinical_tnbc %in% TRUE)
  basal <- dat %>% filter(basal_like_proxy %in% TRUE)
  model_tbl <- bind_rows(
    fit_endpoint(tnbc, "clinical_TNBC", "RFS", "rfs_time_months", "rfs_event"),
    fit_endpoint(tnbc, "clinical_TNBC", "OS", "os_time_months", "os_event"),
    fit_endpoint(basal, "basal_like_proxy", "RFS", "rfs_time_months", "rfs_event"),
    fit_endpoint(basal, "basal_like_proxy", "OS", "os_time_months", "os_event")
  )
  readr::write_tsv(model_tbl, file.path(root, "results", "METABRIC_survival_models.tsv"))

  plot_km(tnbc, "clinical_TNBC", "RFS", "rfs_time_months", "rfs_event",
          file.path(root, "figures", "METABRIC_clinical_TNBC_RFS_BIR_KM.pdf"))
  plot_km(tnbc, "clinical_TNBC", "OS", "os_time_months", "os_event",
          file.path(root, "figures", "METABRIC_clinical_TNBC_OS_BIR_KM.pdf"))
  plot_km(basal, "basal_like_proxy", "RFS", "rfs_time_months", "rfs_event",
          file.path(root, "figures", "METABRIC_basal_like_proxy_RFS_BIR_KM.pdf"))

  message("METABRIC survival validation complete.")
}

main()
