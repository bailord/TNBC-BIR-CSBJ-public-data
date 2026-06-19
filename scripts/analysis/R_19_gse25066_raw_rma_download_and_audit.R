#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(affy)
  library(Biobase)
  library(broom)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(survival)
  library(tibble)
  library(tidyr)
  library(yaml)
})

source("scripts/R/utils_signature_scoring.R")
set.seed(42)

root <- getwd()
cfg <- yaml::read_yaml(file.path(root, "config", "project_config.yaml"))
public_dir <- file.path(cfg$paths$public_data_root, "GEO", "GSE25066")
cel_dir <- file.path(public_dir, "RMA_CEL_cache")
rma_dir <- file.path(public_dir, "RMA_cache")
results_dir <- file.path(root, cfg$paths$results)
fig_dir <- file.path(root, cfg$paths$figures)
reports_dir <- file.path(root, "reports")
logs_dir <- file.path(root, cfg$paths$logs)
processed_dir <- file.path(root, cfg$paths$processed_data)
for (d in c(public_dir, cel_dir, rma_dir, results_dir, fig_dir, reports_dir, logs_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

sig_tbl <- readr::read_tsv(file.path(processed_dir, "TNBC_BIR_signatures_clean.tsv"), show_col_types = FALSE)

raw_manifest <- tibble(
  dataset = c("GSE25055", "GSE25065"),
  url = c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25055/suppl/GSE25055_RAW.tar",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25065/suppl/GSE25065_RAW.tar"
  ),
  public_path = file.path(public_dir, c("GSE25055_RAW.tar", "GSE25065_RAW.tar")),
  expected_bytes = c(770519040, 433418240)
)

sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  out <- tryCatch(system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE), error = function(e) NA_character_)
  if (length(out) == 0 || is.na(out[1])) return(NA_character_)
  strsplit(out[1], "\\s+")[[1]][1]
}

download_one <- function(url, dest) {
  if (file.exists(dest) && file.info(dest)$size > 0) {
    return(tibble(url = url, dest = dest, status = "present", command = NA_character_, message = "Existing non-empty file reused."))
  }
  cmd <- c("--resolve", "ftp.ncbi.nlm.nih.gov:443:130.14.250.13", "-L", "--fail", "--retry", "3", "-C", "-", "-o", dest, url)
  out <- tryCatch(system2("curl", cmd, stdout = TRUE, stderr = TRUE), error = function(e) conditionMessage(e))
  ok <- file.exists(dest) && file.info(dest)$size > 0
  tibble(
    url = url,
    dest = dest,
    status = if (ok) "downloaded" else "failed",
    command = paste("curl", paste(shQuote(cmd), collapse = " ")),
    message = paste(out, collapse = "\n")
  )
}

write_download_blocker <- function(qc) {
  lines <- c(
    "# GSE25066 RAW/RMA Download Blocker",
    "",
    paste("Date:", Sys.Date()),
    "",
    "RAW/RMA could not proceed because one or more RAW tar files failed to download.",
    "",
    "## Attempts",
    "",
    paste0("- ", basename(qc$dest), ": status=", qc$status, "; url=", qc$url),
    "",
    "## Commands",
    "",
    "```sh",
    qc$command[!is.na(qc$command)],
    "```"
  )
  writeLines(lines, file.path(logs_dir, "GSE25066_RMA_download_blocker.md"))
}

extract_tar_if_needed <- function(tar_path, dataset) {
  out_dir <- file.path(cel_dir, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cel_files <- list.files(out_dir, pattern = "[.]CEL([.]gz)?$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(cel_files) > 0) return(cel_files)
  status <- system2("tar", c("-xf", tar_path, "-C", out_dir))
  if (status != 0) stop("tar extraction failed for ", tar_path)
  list.files(out_dir, pattern = "[.]CEL([.]gz)?$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
}

read_gpl_annotation <- function(path) {
  if (!file.exists(path)) stop("Missing GPL annotation: ", path)
  lines <- readLines(gzfile(path), warn = FALSE)
  begin <- grep("^!platform_table_begin", lines)
  end <- grep("^!platform_table_end", lines)
  if (length(begin) == 1 && length(end) == 1) {
    ann <- read.delim(text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"), sep = "\t", header = TRUE, quote = "", comment.char = "", fill = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    ann <- read.delim(gzfile(path), sep = "\t", header = TRUE, comment.char = "#", quote = "", fill = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  }
  names_l <- tolower(names(ann))
  id_col <- names(ann)[match("id", names_l)]
  if (is.na(id_col)) id_col <- names(ann)[grepl("^id$|id_ref|probe", names_l)][1]
  gene_col <- names(ann)[grepl("gene.?symbol|gene symbol|symbol", names_l)][1]
  if (is.na(id_col) || is.na(gene_col)) stop("Could not identify probe/gene columns in GPL annotation.")
  ann %>%
    transmute(
      probe_id = as.character(.data[[id_col]]),
      gene_symbol = toupper(trimws(str_split_fixed(as.character(.data[[gene_col]]), " /// |//|;|,", 2)[, 1]))
    ) %>%
    filter(!is.na(probe_id), probe_id != "", !is.na(gene_symbol), gene_symbol != "", gene_symbol != "---") %>%
    distinct(probe_id, gene_symbol)
}

collapse_probe_matrix <- function(expr_probe, ann, rule = c("highest_iqr", "median")) {
  rule <- match.arg(rule)
  common <- intersect(rownames(expr_probe), ann$probe_id)
  expr_probe <- expr_probe[common, , drop = FALSE]
  ann2 <- ann %>% filter(probe_id %in% common) %>% distinct(probe_id, gene_symbol)
  expr_probe <- expr_probe[ann2$probe_id, , drop = FALSE]
  rownames(expr_probe) <- ann2$gene_symbol
  if (rule == "median") {
    df <- as.data.frame(expr_probe, check.names = FALSE) %>% rownames_to_column("gene")
    out <- df %>% group_by(gene) %>% summarise(across(where(is.numeric), median, na.rm = TRUE), .groups = "drop")
    mat <- as.matrix(out[, -1, drop = FALSE])
    rownames(mat) <- out$gene
    return(mat)
  }
  make_unique_gene_matrix(expr_probe, duplicate_rule = "highest_mad")
}

add_apc_axis <- function(dat) {
  dat %>%
    mutate(
      APC_CXCL9_axis = rowMeans(cbind(
        as.numeric(scale(HLAII)),
        as.numeric(scale(CXCL9_FOLR2_macrophage)),
        as.numeric(scale(DC_HLA_antigen_presentation))
      ), na.rm = TRUE),
      .after = TNBC_BIR_favorable_score
    )
}

fit_logistic <- function(dat, subset_name, model_name, formula) {
  vars <- all.vars(formula)
  md <- dat %>% select(any_of(unique(c("sample_id", vars)))) %>% filter(if_all(all_of(vars), ~ !is.na(.))) %>% mutate(cohort = factor(cohort))
  n_event <- sum(md$pCR_binary == 1)
  n_rd <- sum(md$pCR_binary == 0)
  if (nrow(md) < 20 || n_event < 5 || n_rd < 5) {
    return(tibble(dataset = "GSE25066_RMA", analysis_subset = subset_name, endpoint = "pCR", n_total = nrow(md), n_event_or_pCR = n_event, model_name = model_name, variable = NA_character_, effect_type = "OR", estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_, covariates = paste(deparse(formula), collapse = " "), notes = "Insufficient complete cases/events."))
  }
  fit <- glm(formula, family = binomial(), data = md)
  broom::tidy(fit) %>%
    mutate(beta = estimate) %>%
    transmute(dataset = "GSE25066_RMA", analysis_subset = subset_name, endpoint = "pCR", n_total = nrow(md), n_event_or_pCR = n_event, model_name = model_name, variable = term, effect_type = "OR", estimate = exp(beta), conf_low = exp(beta - 1.96 * std.error), conf_high = exp(beta + 1.96 * std.error), p_value = p.value, covariates = paste(deparse(formula), collapse = " "), notes = "")
}

fit_cox <- function(dat, subset_name, model_name, formula) {
  vars <- all.vars(formula)
  md <- dat %>% select(any_of(unique(c("sample_id", vars)))) %>% filter(if_all(all_of(vars), ~ !is.na(.))) %>% mutate(cohort = factor(cohort))
  n_event <- sum(md$drfs_event == 1)
  if (nrow(md) < 30 || n_event < 10) {
    return(tibble(dataset = "GSE25066_RMA", analysis_subset = subset_name, endpoint = "DRFS", n_total = nrow(md), n_event_or_pCR = n_event, model_name = model_name, variable = NA_character_, effect_type = "HR", estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_, covariates = paste(deparse(formula), collapse = " "), notes = "Insufficient complete cases/events."))
  }
  fit <- coxph(formula, data = md)
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    transmute(dataset = "GSE25066_RMA", analysis_subset = subset_name, endpoint = "DRFS", n_total = nrow(md), n_event_or_pCR = n_event, model_name = model_name, variable = term, effect_type = "HR", estimate = estimate, conf_low = conf.low, conf_high = conf.high, p_value = p.value, covariates = paste(deparse(formula), collapse = " "), notes = "")
}

main <- function() {
  download_qc <- bind_rows(lapply(seq_len(nrow(raw_manifest)), function(i) download_one(raw_manifest$url[i], raw_manifest$public_path[i]))) %>%
    left_join(raw_manifest, by = c("url", "dest" = "public_path")) %>%
    mutate(file_exists = file.exists(dest), file_size_bytes = if_else(file_exists, as.numeric(file.info(dest)$size), NA_real_), sha256 = vapply(dest, sha256_file, character(1)))
  readr::write_tsv(download_qc, file.path(results_dir, "GSE25066_RMA_download_qc.tsv"))
  if (!all(download_qc$file_exists)) {
    write_download_blocker(download_qc)
    stop("GSE25066 RAW download failed; see logs/GSE25066_RMA_download_blocker.md", call. = FALSE)
  }

  cel_files <- unlist(mapply(extract_tar_if_needed, raw_manifest$public_path, raw_manifest$dataset, SIMPLIFY = FALSE), use.names = FALSE)
  cel_qc <- tibble(cel_file = cel_files, sample_id = str_remove(basename(cel_files), "[.](CEL|cel)([.]gz)?$")) %>%
    mutate(sample_id = str_replace(sample_id, "_.*$", ""))
  rma_probe_path <- file.path(rma_dir, "GSE25066_RMA_probe_expression.rds")
  if (file.exists(rma_probe_path)) {
    expr_probe <- readRDS(rma_probe_path)
  } else {
    message("Running affy RMA on ", length(cel_files), " CEL files.")
    affy_obj <- affy::ReadAffy(filenames = cel_files)
    eset <- affy::rma(affy_obj)
    expr_probe <- Biobase::exprs(eset)
    colnames(expr_probe) <- str_remove(basename(sampleNames(eset)), "[.](CEL|cel)([.]gz)?$")
    colnames(expr_probe) <- str_replace(colnames(expr_probe), "_.*$", "")
    saveRDS(expr_probe, rma_probe_path)
  }

  ann <- read_gpl_annotation(file.path(public_dir, "GPL96.annot.gz"))
  expr_gene <- collapse_probe_matrix(expr_probe, ann, "highest_iqr")
  expr_gene_median <- collapse_probe_matrix(expr_probe, ann, "median")
  saveRDS(expr_gene, file.path(rma_dir, "GSE25066_RMA_gene_expression_highest_iqr.rds"))
  saveRDS(expr_gene_median, file.path(rma_dir, "GSE25066_RMA_gene_expression_median.rds"))

  clin <- readr::read_tsv(file.path(results_dir, "GSE25066_TNBC_scores_clinical.tsv"), show_col_types = FALSE) %>%
    select(sample_id, geo_accession, patient_title, cohort, ER_raw, PR_raw, HER2_raw, pCR_raw, DRFS_event_raw, DRFS_time_raw, PAM50_raw, age, clinical_stage, grade, ER_binary, PR_binary, HER2_binary, ER_negative, PR_negative, HER2_negative, primary_tnbc, ER_HER2_negative_sensitivity, basal_like_sensitivity, pCR_binary, drfs_event, drfs_time_years, tnbc_definition, endpoint_definition)
  common <- intersect(colnames(expr_gene), clin$sample_id)
  expr_gene <- expr_gene[, common, drop = FALSE]
  clin <- clin %>% filter(sample_id %in% common)

  score_rma <- score_signatures(expr_gene, sig_tbl, min_coverage = cfg$score_settings$min_signature_coverage, duplicate_rule = cfg$score_settings$duplicate_gene_rule)
  scores <- assign_bir_quadrants(score_rma$scores) %>% add_apc_axis()
  scores_clin <- scores %>% left_join(clin, by = "sample_id")
  readr::write_tsv(scores_clin, file.path(results_dir, "GSE25066_RMA_scores_clinical.tsv"))
  readr::write_tsv(score_rma$coverage %>% mutate(dataset = "GSE25066_RMA", .before = 1), file.path(results_dir, "GSE25066_RMA_signature_coverage.tsv"))

  processed <- readr::read_tsv(file.path(results_dir, "GSE25066_TNBC_scores_clinical.tsv"), show_col_types = FALSE) %>% add_apc_axis()
  score_vars <- c("TNBC_BIR_favorable_score", "Immune_remodeling_score", "HLAII", "Proliferation", "APC_CXCL9_axis")
  corr <- lapply(score_vars, function(v) {
    joined <- processed %>% select(sample_id, processed = all_of(v)) %>% inner_join(scores_clin %>% select(sample_id, RMA = all_of(v)), by = "sample_id")
    tibble(score = v, n = nrow(joined), pearson = cor(joined$processed, joined$RMA, use = "complete.obs", method = "pearson"), spearman = cor(joined$processed, joined$RMA, use = "complete.obs", method = "spearman"))
  }) %>% bind_rows() %>% mutate(decision = if_else(pearson >= 0.85 & spearman >= 0.85, "high_correlation", "low_or_mixed_correlation"))
  readr::write_tsv(corr, file.path(results_dir, "GSE25066_processed_vs_RMA_score_correlation.tsv"))

  plot_dat <- bind_rows(lapply(score_vars[1:4], function(v) {
    processed %>% select(sample_id, processed = all_of(v)) %>% inner_join(scores_clin %>% select(sample_id, RMA = all_of(v)), by = "sample_id") %>% mutate(score = v)
  }))
  p_scatter <- ggplot(plot_dat, aes(x = processed, y = RMA)) +
    geom_point(alpha = 0.65, size = 0.9) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = "#B4443E") +
    facet_wrap(~score, scales = "free") +
    theme_bw() +
    labs(x = "Processed-series score", y = "RAW/RMA score", title = "GSE25066 processed vs RAW/RMA score concordance")
  ggsave(file.path(fig_dir, "GSE25066_processed_vs_RMA_score_scatter.pdf"), p_scatter, width = 8.5, height = 6.2)

  primary <- scores_clin %>% filter(primary_tnbc %in% TRUE)
  models <- bind_rows(
    fit_logistic(primary, "primary_TNBC", "pCR_BIR", pCR_binary ~ cohort + scale(TNBC_BIR_favorable_score)),
    fit_logistic(primary, "primary_TNBC", "pCR_Immune", pCR_binary ~ cohort + scale(Immune_remodeling_score)),
    fit_logistic(primary, "primary_TNBC", "pCR_APC_CXCL9_axis", pCR_binary ~ cohort + scale(APC_CXCL9_axis)),
    fit_logistic(primary, "primary_TNBC", "pCR_APC_CXCL9_plus_proliferation", pCR_binary ~ cohort + scale(APC_CXCL9_axis) + scale(Proliferation)),
    fit_logistic(primary, "primary_TNBC", "pCR_Proliferation", pCR_binary ~ cohort + scale(Proliferation)),
    fit_cox(primary, "primary_TNBC", "DRFS_APC_CXCL9_axis", Surv(drfs_time_years, drfs_event) ~ cohort + scale(APC_CXCL9_axis)),
    fit_cox(primary, "primary_TNBC", "DRFS_Immune", Surv(drfs_time_years, drfs_event) ~ cohort + scale(Immune_remodeling_score))
  )
  readr::write_tsv(models, file.path(results_dir, "GSE25066_RMA_pCR_DRFS_models.tsv"))

  qc <- tibble(
    dataset = "GSE25066",
    metric = c("raw_tar_files", "cel_files", "rma_probe_features", "rma_samples", "gene_features_highest_iqr", "clinical_join_samples", "primary_tnbc_samples", "primary_tnbc_pCR", "primary_tnbc_DRFS_events"),
    value = c(nrow(raw_manifest), length(cel_files), nrow(expr_probe), ncol(expr_probe), nrow(expr_gene), nrow(scores_clin), nrow(primary), sum(primary$pCR_binary == 1, na.rm = TRUE), sum(primary$drfs_event == 1, na.rm = TRUE)),
    notes = c("Expected GSE25055/GSE25065 RAW tar files.", "CEL files extracted from RAW tar cache.", "Probe rows after affy RMA.", "RMA sample columns.", "Gene rows after GPL96 mapping and highest-IQR collapse.", "Samples joined to processed clinical mapping.", "Primary ER-/PR-/HER2- subset.", "pCR events in primary TNBC.", "DRFS events in primary TNBC.")
  )
  readr::write_tsv(qc, file.path(results_dir, "GSE25066_RMA_expression_qc.tsv"))

  unchanged <- all(corr$pearson >= 0.85 & corr$spearman >= 0.85, na.rm = TRUE)
  lines <- c(
    "# GSE25066 RAW/RMA Audit Report",
    "",
    paste("Generated:", Sys.time()),
    "",
    paste0("RAW/RMA status: ", if (unchanged) "processed-series scores are technically concordant with RAW/RMA." else "processed-series and RAW/RMA scores are discordant or mixed."),
    "",
    "## Score Concordance",
    "",
    paste0("- ", corr$score, ": Pearson=", sprintf("%.3f", corr$pearson), ", Spearman=", sprintf("%.3f", corr$spearman), " (", corr$decision, ")"),
    "",
    "## Decision",
    "",
    if (unchanged) "Accept processed GSE25066 as technically adequate for the negative external pCR conclusion." else "Use caution: inspect RMA models and prefer RMA if conclusions differ from processed results.",
    "",
    "## Files",
    "",
    "- `results/GSE25066_processed_vs_RMA_score_correlation.tsv`",
    "- `results/GSE25066_RMA_pCR_DRFS_models.tsv`",
    "- `figures/GSE25066_processed_vs_RMA_score_scatter.pdf`"
  )
  writeLines(lines, file.path(reports_dir, "GSE25066_RMA_audit_report.md"))
  message("GSE25066 RAW/RMA audit complete.")
}

main()
