#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(ggplot2)
  library(purrr)
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
raw_dir <- file.path(root, cfg$paths$raw_data, "GSE25066")
processed_dir <- file.path(root, cfg$paths$processed_data)
results_dir <- file.path(root, cfg$paths$results)
fig_dir <- file.path(root, cfg$paths$figures)
log_dir <- file.path(root, cfg$paths$logs)
for (d in c(public_dir, raw_dir, processed_dir, results_dir, fig_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

sig_tbl <- readr::read_tsv(file.path(processed_dir, "TNBC_BIR_signatures_clean.tsv"), show_col_types = FALSE)

sha256_file <- function(path) {
  out <- tryCatch(system2("shasum", c("-a", "256", path), stdout = TRUE), error = function(e) NA_character_)
  if (length(out) == 0 || is.na(out[1])) return(NA_character_)
  strsplit(out[1], "\\s+")[[1]][1]
}

append_manifest <- function(paths, phase = "Phase4_GSE25066_external_validation") {
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) return(invisible(FALSE))
  manifest_path <- file.path(results_dir, "public_data_acquisition_manifest.tsv")
  rows <- tibble(
    dataset = "GSE25066",
    source_url = "NCBI GEO processed series matrix or GPL96 annotation, reused from public data root when present.",
    public_path = normalizePath(paths, mustWork = TRUE),
    sha256 = vapply(paths, sha256_file, character(1)),
    download_date = format(Sys.Date()),
    used_in_phase = phase
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
  invisible(TRUE)
}

download_if_missing <- function(urls, dest) {
  if (file.exists(dest) && file.info(dest)$size > 0) return(TRUE)
  for (url in urls) {
    ok <- tryCatch({
      utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) FALSE)
    if (!ok && nzchar(Sys.which("curl"))) {
      status <- suppressWarnings(system2("curl", c("-L", "--fail", "--retry", "2", "-o", dest, url)))
      ok <- status == 0
    }
    if (!ok && nzchar(Sys.which("curl")) && grepl("https://ftp[.]ncbi[.]nlm[.]nih[.]gov", url)) {
      status <- suppressWarnings(system2(
        "curl",
        c(
          "--resolve", "ftp.ncbi.nlm.nih.gov:443:130.14.250.13",
          "-L", "--fail", "--retry", "2", "-C", "-", "-o", dest, url
        )
      ))
      ok <- status == 0
    }
    if (ok && file.exists(dest) && file.info(dest)$size > 0) return(TRUE)
    if (file.exists(dest)) file.remove(dest)
  }
  FALSE
}

link_into_project <- function(public_path) {
  project_path <- file.path(raw_dir, basename(public_path))
  if (file.exists(project_path) || nzchar(Sys.readlink(project_path))) file.remove(project_path)
  ok <- file.symlink(public_path, project_path)
  if (!ok) file.copy(public_path, project_path, overwrite = TRUE)
  invisible(project_path)
}

write_blocker_outputs <- function(reason) {
  blocker <- c(
    "# GSE25066 Data Acquisition / Mapping Blocker",
    "",
    paste("Date:", Sys.Date()),
    "",
    "## Status",
    "",
    reason,
    "",
    "## Required public files",
    "",
    "- `GSE25055_series_matrix.txt.gz`",
    "- `GSE25065_series_matrix.txt.gz`",
    "- GPL96 annotation, preferably `GPL96.annot.gz` or `GPL96_family.soft.gz`",
    "",
    "Place these files under `<PUBLIC_BIOINF_DATABASES>/GEO/GSE25066/` and rerun:",
    "",
    "```sh",
    "Rscript scripts/R/10_gse25066_external_validation.R",
    "```"
  )
  writeLines(blocker, file.path(log_dir, "GSE25066_data_acquisition_blocker.md"))

  readr::write_tsv(tibble(
    dataset = "GSE25066",
    step = c("processed_matrix_available", "gpl96_annotation_available", "analysis_status"),
    n = c(0, 0, 0),
    notes = c("GSE25055/GSE25065 processed matrix files are not available locally and NCBI download failed.",
              "GPL96 annotation is not available locally and NCBI download failed.",
              reason)
  ), file.path(results_dir, "GSE25066_sample_count_flow.tsv"))
  readr::write_tsv(tibble(
    dataset = "GSE25066",
    field = c("pathologic_response_pcr_rd", "drfs_1_event_0_censored", "drfs_even_time_years"),
    mapped_name = c("pCR_binary", "drfs_event", "drfs_time_years"),
    n_missing = NA_integer_,
    missing_percent = NA_real_,
    status = "not_mapped",
    notes = reason
  ), file.path(results_dir, "GSE25066_endpoint_mapping_qc.tsv"))
  readr::write_tsv(tibble(dataset = "GSE25066", signature = NA_character_, n_signature_genes = NA_integer_,
                          n_genes_found = NA_integer_, coverage = NA_real_, missing_genes = NA_character_,
                          notes = reason),
                   file.path(results_dir, "GSE25066_signature_coverage.tsv"))
  readr::write_tsv(tibble(), file.path(results_dir, "GSE25066_TNBC_scores_clinical.tsv"))
  readr::write_tsv(tibble(dataset = "GSE25066", model_name = NA_character_, notes = reason),
                   file.path(results_dir, "GSE25066_pCR_models.tsv"))
  readr::write_tsv(tibble(dataset = "GSE25066", model_name = NA_character_, notes = reason),
                   file.path(results_dir, "GSE25066_pCR_model_auc.tsv"))
  readr::write_tsv(tibble(dataset = "GSE25066", comparison = NA_character_, notes = reason),
                   file.path(results_dir, "GSE25066_pCR_model_lrt.tsv"))
  readr::write_tsv(tibble(dataset = "GSE25066", model_name = NA_character_, notes = reason),
                   file.path(results_dir, "GSE25066_DRFS_models.tsv"))
  readr::write_tsv(tibble(dataset = "GSE25066", model_name = NA_character_, notes = reason),
                   file.path(results_dir, "GSE25066_DRFS_PH_assumption.tsv"))
  stop(reason, call. = FALSE)
}

ensure_public_files <- function() {
  matrix_files <- file.path(public_dir, c("GSE25055_series_matrix.txt.gz", "GSE25065_series_matrix.txt.gz"))
  gpl_files <- file.path(public_dir, c("GPL96.annot.gz", "GPL96_family.soft.gz"))

  base <- "https://ftp.ncbi.nlm.nih.gov/geo"
  urls <- list(
    GSE25055 = c(
      paste0(base, "/series/GSE25nnn/GSE25055/matrix/GSE25055_series_matrix.txt.gz"),
      "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE25055&format=file&file=GSE25055_series_matrix.txt.gz"
    ),
    GSE25065 = c(
      paste0(base, "/series/GSE25nnn/GSE25065/matrix/GSE25065_series_matrix.txt.gz"),
      "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE25065&format=file&file=GSE25065_series_matrix.txt.gz"
    ),
    GPL96_annot = c(
      paste0(base, "/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz"),
      "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GPL96&format=file&file=GPL96.annot.gz"
    ),
    GPL96_soft = c(
      paste0(base, "/platforms/GPLnnn/GPL96/soft/GPL96_family.soft.gz"),
      "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GPL96&format=file&file=GPL96_family.soft.gz"
    )
  )
  download_if_missing(urls$GSE25055, matrix_files[1])
  download_if_missing(urls$GSE25065, matrix_files[2])
  if (!any(file.exists(gpl_files) & file.info(gpl_files)$size > 0)) {
    download_if_missing(urls$GPL96_annot, gpl_files[1])
    if (!file.exists(gpl_files[1])) download_if_missing(urls$GPL96_soft, gpl_files[2])
  }

  matrix_ok <- file.exists(matrix_files) & file.info(matrix_files)$size > 0
  gpl_ok <- file.exists(gpl_files) & file.info(gpl_files)$size > 0
  if (!all(matrix_ok) || !any(gpl_ok)) {
    write_blocker_outputs("GSE25066 processed GEO files are not available locally and this machine cannot download them from NCBI/GEO due TLS/connection failure.")
  }
  lapply(c(matrix_files, gpl_files[gpl_ok]), link_into_project)
  append_manifest(c(matrix_files, gpl_files[gpl_ok]))
  list(matrix_files = matrix_files, gpl_file = gpl_files[gpl_ok][1])
}

parse_matrix_line <- function(line) {
  read.delim(text = line, sep = "\t", header = FALSE, quote = "\"", comment.char = "",
             check.names = FALSE, stringsAsFactors = FALSE)
}

parse_series_matrix <- function(path, cohort) {
  message("Parsing ", cohort, ": ", path)
  lines <- readLines(gzfile(path), warn = FALSE)
  sample_lines <- lines[grepl("^!Sample_", lines)]
  parsed <- lapply(sample_lines, parse_matrix_line)
  get_values <- function(key) {
    hits <- parsed[vapply(parsed, function(x) identical(x[1, 1], key), logical(1))]
    if (length(hits) == 0) return(NULL)
    as.character(hits[[1]][1, -1])
  }
  geo <- get_values("!Sample_geo_accession")
  title <- get_values("!Sample_title")
  if (is.null(geo)) stop("Could not find !Sample_geo_accession in ", path)
  pheno <- tibble(
    sample_id = geo,
    geo_accession = geo,
    patient_title = if (!is.null(title)) title else geo,
    cohort = cohort
  )
  characteristic_rows <- parsed[vapply(parsed, function(x) identical(x[1, 1], "!Sample_characteristics_ch1"), logical(1))]
  for (row in characteristic_rows) {
    vals <- as.character(row[1, -1])
    keys <- sub(":.*$", "", vals)
    key <- names(sort(table(keys), decreasing = TRUE))[1]
    col <- key %>%
      str_replace_all("[^A-Za-z0-9]+", "_") %>%
      str_replace_all("^_|_$", "") %>%
      tolower()
    value <- ifelse(grepl(":", vals), sub("^[^:]+:\\s*", "", vals), vals)
    if (col %in% names(pheno)) col <- make.unique(c(names(pheno), col))[length(names(pheno)) + 1]
    pheno[[col]] <- value
  }

  begin <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end", lines)
  if (length(begin) != 1 || length(end) != 1) stop("Could not find series matrix table boundaries in ", path)
  tbl <- read.delim(text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"),
                    sep = "\t", header = TRUE, quote = "\"", comment.char = "",
                    check.names = FALSE, stringsAsFactors = FALSE)
  probe <- tbl[[1]]
  expr <- as.matrix(as.data.frame(lapply(tbl[, -1, drop = FALSE], function(x) suppressWarnings(as.numeric(x))), check.names = FALSE))
  rownames(expr) <- probe
  colnames(expr) <- colnames(tbl)[-1]
  common <- intersect(colnames(expr), pheno$sample_id)
  expr <- expr[, common, drop = FALSE]
  pheno <- pheno %>% filter(sample_id %in% common)
  list(expr = expr, pheno = pheno)
}

read_gpl96_annotation <- function(path) {
  message("Reading GPL96 annotation: ", path)
  lines <- readLines(gzfile(path), warn = FALSE)
  begin <- grep("^!platform_table_begin", lines)
  end <- grep("^!platform_table_end", lines)
  if (length(begin) == 1 && length(end) == 1) {
    ann <- read.delim(text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"),
                      sep = "\t", header = TRUE, quote = "", comment.char = "",
                      fill = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (grepl("annot[.]gz$", path, ignore.case = TRUE)) {
    ann <- read.delim(gzfile(path), sep = "\t", header = TRUE, comment.char = "#",
                      quote = "", fill = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    stop("Could not find GPL96 platform table boundaries.")
  }
  names_l <- tolower(names(ann))
  id_col <- names(ann)[match("id", names_l)]
  if (is.na(id_col)) id_col <- names(ann)[grepl("^id$|id_ref|probe", names_l)][1]
  gene_col <- names(ann)[grepl("gene.?symbol|gene symbol|symbol", names_l)][1]
  if (is.na(id_col) || is.na(gene_col)) {
    stop("Could not identify probe ID and gene symbol columns in GPL96 annotation. Columns: ", paste(names(ann), collapse = ", "))
  }
  ann %>%
    transmute(
      probe_id = as.character(.data[[id_col]]),
      gene_symbol = toupper(trimws(str_split_fixed(as.character(.data[[gene_col]]), " /// |//|;|,", 2)[, 1]))
    ) %>%
    filter(!is.na(probe_id), probe_id != "", !is.na(gene_symbol), gene_symbol != "", gene_symbol != "---") %>%
    distinct(probe_id, gene_symbol)
}

as_binary_status <- function(x, positive_words = c("positive", "pos", "1"), negative_words = c("negative", "neg", "0")) {
  y <- tolower(str_squish(as.character(x)))
  case_when(
    y %in% c(positive_words, "p", "pos", "+") | str_detect(y, "positive") ~ 1,
    y %in% c(negative_words, "n", "neg", "-") | str_detect(y, "negative") ~ 0,
    y %in% c("i", "indeterminate", "equivocal", "unknown", "na", "n/a", "") ~ NA_real_,
    TRUE ~ NA_real_
  )
}

map_pcr <- function(x) {
  y <- toupper(str_squish(as.character(x)))
  case_when(
    y == "PCR" | str_detect(y, "^PCR$|P\\.CR|COMPLETE") ~ 1,
    y == "RD" | str_detect(y, "RESIDUAL|NON") ~ 0,
    TRUE ~ NA_real_
  )
}

find_col <- function(dat, candidates) {
  nms <- names(dat)
  nms_l <- tolower(nms)
  cand_l <- tolower(candidates)
  hit <- nms[match(cand_l, nms_l, nomatch = 0)]
  hit <- hit[hit != ""]
  if (length(hit) > 0) hit[1] else NA_character_
}

map_clinical <- function(pheno) {
  er_col <- find_col(pheno, c("er_status_ihc", "esr1_status", "er_status"))
  pr_col <- find_col(pheno, c("pr_status_ihc", "pr_status"))
  her2_col <- find_col(pheno, c("her2_status", "erbb2_status"))
  pcr_col <- find_col(pheno, c("pathologic_response_pcr_rd", "pcr_rd", "pathologic_response"))
  drfs_event_col <- find_col(pheno, c("drfs_1_event_0_censored", "drfs_event"))
  drfs_time_col <- find_col(pheno, c("drfs_even_time_years", "drfs_event_time_years", "drfs_time_years"))
  pam50_col <- find_col(pheno, c("pam50_class", "pam50", "pam50_subtype"))
  age_col <- find_col(pheno, c("age", "age_at_diagnosis"))
  stage_col <- find_col(pheno, c("clinical_stage", "stage", "tumor_stage"))
  grade_col <- find_col(pheno, c("grade", "tumor_grade"))

  mapped <- pheno %>%
    transmute(
      sample_id, geo_accession, patient_title, cohort,
      ER_raw = if (!is.na(er_col)) .data[[er_col]] else NA_character_,
      PR_raw = if (!is.na(pr_col)) .data[[pr_col]] else NA_character_,
      HER2_raw = if (!is.na(her2_col)) .data[[her2_col]] else NA_character_,
      pCR_raw = if (!is.na(pcr_col)) .data[[pcr_col]] else NA_character_,
      DRFS_event_raw = if (!is.na(drfs_event_col)) .data[[drfs_event_col]] else NA_character_,
      DRFS_time_raw = if (!is.na(drfs_time_col)) .data[[drfs_time_col]] else NA_character_,
      PAM50_raw = if (!is.na(pam50_col)) .data[[pam50_col]] else NA_character_,
      age = if (!is.na(age_col)) suppressWarnings(as.numeric(.data[[age_col]])) else NA_real_,
      clinical_stage = if (!is.na(stage_col)) as.character(.data[[stage_col]]) else NA_character_,
      grade = if (!is.na(grade_col)) as.character(.data[[grade_col]]) else NA_character_
    ) %>%
    mutate(
      ER_binary = as_binary_status(ER_raw),
      PR_binary = as_binary_status(PR_raw),
      HER2_binary = as_binary_status(HER2_raw),
      ER_negative = ER_binary == 0,
      PR_negative = PR_binary == 0,
      HER2_negative = HER2_binary == 0,
      primary_tnbc = ER_negative & PR_negative & HER2_negative,
      ER_HER2_negative_sensitivity = ER_negative & HER2_negative,
      basal_like_sensitivity = str_detect(PAM50_raw, regex("^Basal", ignore_case = TRUE)),
      pCR_binary = map_pcr(pCR_raw),
      drfs_event = suppressWarnings(as.numeric(DRFS_event_raw)),
      drfs_time_years = suppressWarnings(as.numeric(DRFS_time_raw)),
      tnbc_definition = "Primary: ER-/PR-/HER2- from GEO IHC/status fields; sensitivity: ER-/HER2- and PAM50 basal-like.",
      endpoint_definition = "pCR/RD from pathologic_response_pcr_rd; DRFS from drfs_1_event_0_censored and drfs_even_time_years."
    )
  attr(mapped, "mapping_cols") <- tibble(
    endpoint_or_covariate = c("ER", "PR", "HER2", "pCR", "DRFS_event", "DRFS_time", "PAM50", "age", "stage", "grade"),
    raw_field_used = c(er_col, pr_col, her2_col, pcr_col, drfs_event_col, drfs_time_col, pam50_col, age_col, stage_col, grade_col)
  )
  mapped
}

endpoint_qc <- function(clin, mapping_cols) {
  fields <- tibble(
    field = c("pathologic_response_pcr_rd", "drfs_1_event_0_censored", "drfs_even_time_years"),
    mapped_name = c("pCR_binary", "drfs_event", "drfs_time_years")
  )
  fields %>%
    mutate(
      dataset = "GSE25066",
      n_missing = vapply(mapped_name, function(x) sum(is.na(clin[[x]])), integer(1)),
      missing_percent = n_missing / nrow(clin),
      raw_field_used = c(
        mapping_cols$raw_field_used[mapping_cols$endpoint_or_covariate == "pCR"],
        mapping_cols$raw_field_used[mapping_cols$endpoint_or_covariate == "DRFS_event"],
        mapping_cols$raw_field_used[mapping_cols$endpoint_or_covariate == "DRFS_time"]
      ),
      status = ifelse(is.na(raw_field_used) | missing_percent > 0.20, "review_required", "mapped"),
      notes = "Stop for review if pCR or DRFS metadata are unmapped for >20% of processed samples."
    ) %>%
    select(dataset, field, raw_field_used, mapped_name, n_missing, missing_percent, status, notes)
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

fit_pcr_models <- function(dat) {
  model_dat <- dat %>%
    filter(primary_tnbc %in% TRUE, !is.na(pCR_binary), !is.na(TNBC_BIR_favorable_score)) %>%
    mutate(cohort = factor(cohort))
  specs <- list(
    Model0_cohort = pCR_binary ~ cohort,
    Model1_cohort_BIR = pCR_binary ~ cohort + scale(TNBC_BIR_favorable_score),
    Model2_cohort_components = pCR_binary ~ cohort + scale(Basoluminal_score) + scale(Immune_remodeling_score) + scale(Suppressive_myeloid_score),
    Immune_only = pCR_binary ~ cohort + scale(Immune_remodeling_score)
  )
  if (nrow(model_dat) < 20 || sum(model_dat$pCR_binary == 1) < 5 || sum(model_dat$pCR_binary == 0) < 5) {
    note <- "Insufficient primary TNBC pCR samples/events for stable pCR modeling."
    return(list(model_tbl = tibble(dataset = "GSE25066", model_name = NA_character_, notes = note),
                auc_tbl = tibble(dataset = "GSE25066", model_name = NA_character_, notes = note),
                lrt_tbl = tibble(dataset = "GSE25066", comparison = NA_character_, notes = note),
                pred_tbl = tibble()))
  }
  mods <- purrr::map(specs, ~ glm(.x, family = binomial(), data = model_dat))
  model_tbl <- purrr::imap_dfr(mods, function(m, nm) {
    broom::tidy(m) %>%
      mutate(
        dataset = "GSE25066",
        analysis_subset = "primary_TNBC_ER_PR_HER2_negative",
        endpoint = "pCR",
        n_total = nobs(m),
        n_event_or_pCR = sum(model_dat$pCR_binary == 1),
        model_name = nm,
        variable = term,
        effect_type = "OR",
        estimate = exp(estimate),
        conf_low = exp(log(estimate) - 1.96 * std.error),
        conf_high = exp(log(estimate) + 1.96 * std.error),
        p_value = p.value,
        covariates = paste(deparse(specs[[nm]]), collapse = " "),
        notes = "GSE25066 processed matrix primary TNBC pCR model.",
        .keep = "unused"
      ) %>%
      select(dataset, analysis_subset, endpoint, n_total, n_event_or_pCR, model_name, variable,
             effect_type, estimate, conf_low, conf_high, p_value, covariates, notes)
  })
  pred_tbl <- purrr::imap_dfr(mods, function(m, nm) {
    tibble(sample_id = model_dat$sample_id, model_name = nm, pCR_binary = model_dat$pCR_binary, pred = as.numeric(predict(m, type = "response")))
  })
  auc_tbl <- pred_tbl %>%
    group_by(model_name) %>%
    summarise(dataset = "GSE25066", auc = calc_auc(pCR_binary, pred), n = n(),
              n_pCR = sum(pCR_binary == 1), n_residual = sum(pCR_binary == 0), .groups = "drop") %>%
    select(dataset, model_name, auc, n, n_pCR, n_residual)
  lrt_tbl <- tibble(
    dataset = "GSE25066",
    comparison = c("Model1_vs_Model0", "Model2_vs_Model0", "Immune_only_vs_Model0"),
    p_value_lrt = c(
      anova(mods$Model0_cohort, mods$Model1_cohort_BIR, test = "Chisq")$`Pr(>Chi)`[2],
      anova(mods$Model0_cohort, mods$Model2_cohort_components, test = "Chisq")$`Pr(>Chi)`[2],
      anova(mods$Model0_cohort, mods$Immune_only, test = "Chisq")$`Pr(>Chi)`[2]
    )
  )
  list(model_tbl = model_tbl, auc_tbl = auc_tbl, lrt_tbl = lrt_tbl, pred_tbl = pred_tbl)
}

fit_drfs_models <- function(dat) {
  model_dat <- dat %>%
    filter(primary_tnbc %in% TRUE, !is.na(drfs_time_years), !is.na(drfs_event), !is.na(TNBC_BIR_favorable_score)) %>%
    mutate(cohort = factor(cohort))
  specs <- list(
    DRFS_BIR = Surv(drfs_time_years, drfs_event) ~ scale(TNBC_BIR_favorable_score) + cohort,
    DRFS_components = Surv(drfs_time_years, drfs_event) ~ scale(Basoluminal_score) + scale(Immune_remodeling_score) + scale(Suppressive_myeloid_score) + cohort,
    DRFS_immune_only = Surv(drfs_time_years, drfs_event) ~ scale(Immune_remodeling_score) + cohort
  )
  if (nrow(model_dat) < 30 || sum(model_dat$drfs_event == 1) < 10) {
    note <- "Insufficient primary TNBC DRFS samples/events for stable Cox modeling."
    return(list(model_tbl = tibble(dataset = "GSE25066", model_name = NA_character_, notes = note),
                ph_tbl = tibble(dataset = "GSE25066", model_name = NA_character_, notes = note)))
  }
  mods <- purrr::map(specs, ~ coxph(.x, data = model_dat))
  model_tbl <- purrr::imap_dfr(mods, function(m, nm) {
    broom::tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
      transmute(
        dataset = "GSE25066",
        analysis_subset = "primary_TNBC_ER_PR_HER2_negative",
        endpoint = "DRFS",
        n_total = nrow(model_dat),
        n_event_or_pCR = sum(model_dat$drfs_event == 1),
        model_name = nm,
        variable = term,
        effect_type = "HR",
        estimate = estimate,
        conf_low = conf.low,
        conf_high = conf.high,
        p_value = p.value,
        covariates = paste(deparse(specs[[nm]]), collapse = " "),
        notes = "Cox model using DRFS years; cohort fixed effect included."
      )
  })
  ph_tbl <- purrr::imap_dfr(mods, function(m, nm) {
    ph <- tryCatch(cox.zph(m), error = function(e) NULL)
    if (is.null(ph)) return(tibble(dataset = "GSE25066", model_name = nm, variable = NA_character_, chisq = NA_real_, p_value = NA_real_, notes = "cox.zph failed."))
    as.data.frame(ph$table) %>%
      rownames_to_column("variable") %>%
      transmute(dataset = "GSE25066", model_name = nm, variable = variable, chisq = chisq, p_value = p, notes = "Schoenfeld residual PH check.")
  })
  list(model_tbl = model_tbl, ph_tbl = ph_tbl, model_data = model_dat)
}

make_figures <- function(dat, pcr_res, drfs_res) {
  tnbc <- dat %>% filter(primary_tnbc %in% TRUE, !is.na(pCR_binary))
  if (nrow(tnbc) > 0) {
    p1 <- tnbc %>%
      mutate(pCR_label = if_else(pCR_binary == 1, "pCR", "Residual disease")) %>%
      ggplot(aes(x = pCR_label, y = Immune_remodeling_score, fill = pCR_label)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.55) +
      geom_jitter(width = 0.15, alpha = 0.65, size = 1.4) +
      theme_bw() +
      theme(legend.position = "none") +
      labs(x = NULL, y = "Immune remodeling score", title = "GSE25066 TNBC: immune remodeling by pCR")
    ggsave(file.path(fig_dir, "GSE25066_scores_by_pCR.pdf"), p1, width = 5.4, height = 4.2)
  }
  if (nrow(pcr_res$pred_tbl) > 0) {
    roc_tbl <- pcr_res$pred_tbl %>%
      group_by(model_name) %>%
      group_modify(~ roc_points(.x$pCR_binary, .x$pred)) %>%
      ungroup() %>%
      left_join(pcr_res$auc_tbl %>% select(model_name, auc), by = "model_name") %>%
      mutate(model_label = paste0(model_name, " (AUC=", sprintf("%.3f", auc), ")"))
    p2 <- roc_tbl %>%
      ggplot(aes(x = 1 - specificity, y = sensitivity, color = model_label)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray55") +
      geom_step(linewidth = 0.9) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      theme_bw() +
      theme(legend.position = "bottom") +
      labs(x = "1 - specificity", y = "Sensitivity", color = NULL, title = "GSE25066 TNBC pCR ROC")
    ggsave(file.path(fig_dir, "GSE25066_ROC_models.pdf"), p2, width = 6.4, height = 5.2)
  }
  if ("model_data" %in% names(drfs_res) && nrow(drfs_res$model_data) > 0) {
    km_dat <- drfs_res$model_data %>%
      mutate(
        immune_tertile = ntile(Immune_remodeling_score, 3),
        immune_group = factor(c("Low immune", "Intermediate", "High immune")[immune_tertile],
                              levels = c("Low immune", "Intermediate", "High immune"))
      )
    fit <- survfit(Surv(drfs_time_years, drfs_event) ~ immune_group, data = km_dat)
    pdf(file.path(fig_dir, "GSE25066_immune_DRFS_KM.pdf"), width = 6.4, height = 5.2)
    plot(fit, col = c("#0072B2", "#666666", "#D55E00"), lwd = 2,
         xlab = "Years", ylab = "DRFS probability", main = "GSE25066 TNBC DRFS by immune remodeling tertile")
    legend("bottomleft", legend = levels(km_dat$immune_group), col = c("#0072B2", "#666666", "#D55E00"), lwd = 2, bty = "n")
    dev.off()
  }
}

main <- function() {
  files <- ensure_public_files()
  gse55 <- parse_series_matrix(files$matrix_files[1], "GSE25055")
  gse65 <- parse_series_matrix(files$matrix_files[2], "GSE25065")
  expr_probe <- cbind(gse55$expr, gse65$expr)
  pheno <- bind_rows(gse55$pheno, gse65$pheno)
  ann <- read_gpl96_annotation(files$gpl_file)

  ann <- ann %>% filter(probe_id %in% rownames(expr_probe))
  expr_gene <- expr_probe[ann$probe_id, , drop = FALSE]
  rownames(expr_gene) <- ann$gene_symbol
  expr_gene <- make_unique_gene_matrix(expr_gene, duplicate_rule = cfg$score_settings$duplicate_gene_rule)

  res <- score_signatures(
    expr_gene,
    sig_tbl = sig_tbl,
    min_coverage = cfg$score_settings$min_signature_coverage,
    duplicate_rule = cfg$score_settings$duplicate_gene_rule
  )
  scores <- assign_bir_quadrants(res$scores)
  coverage <- res$coverage %>% mutate(dataset = "GSE25066", .before = 1)
  min_main_cov <- coverage %>%
    filter(signature %in% c("Basoluminal_plasticity", "Mast_cell", "Tfh_TLS_like", "CD8_Tsem_memory",
                            "CXCL9_FOLR2_macrophage", "DC_HLA_antigen_presentation", "Suppressive_myeloid")) %>%
    summarise(x = min(coverage, na.rm = TRUE)) %>%
    pull(x)
  readr::write_tsv(coverage, file.path(results_dir, "GSE25066_signature_coverage.tsv"))
  if (is.na(min_main_cov) || min_main_cov < 0.70) {
    write_blocker_outputs("GSE25066 processed matrix loaded, but main-module signature coverage is <0.70; do not interpret and consider RAW CEL/RMA.")
  }

  clin <- map_clinical(pheno)
  mapping_cols <- attr(clin, "mapping_cols")
  qc <- endpoint_qc(clin, mapping_cols)
  readr::write_tsv(qc, file.path(results_dir, "GSE25066_endpoint_mapping_qc.tsv"))
  if (any(qc$status == "review_required" & qc$field %in% c("pathologic_response_pcr_rd", "drfs_1_event_0_censored", "drfs_even_time_years"))) {
    readr::write_tsv(scores %>% left_join(clin, by = "sample_id"), file.path(results_dir, "GSE25066_TNBC_scores_clinical.tsv"))
    write_blocker_outputs("GSE25066 endpoint metadata missing/unmapped for >20% of samples; stop for field mapping review.")
  }

  dat <- scores %>% left_join(clin, by = "sample_id")
  readr::write_tsv(dat, file.path(results_dir, "GSE25066_TNBC_scores_clinical.tsv"))
  saveRDS(list(expr_gene = expr_gene, pheno = pheno, scores_clinical = dat, coverage = coverage),
          file.path(processed_dir, "GSE25066_scored_bulk.rds"))

  count_flow <- tibble(
    dataset = "GSE25066",
    step = c("processed_expression_samples", "clinical_samples", "expression_clinical_overlap",
             "primary_TNBC_samples", "ER_HER2_negative_sensitivity_samples", "PAM50_basal_like_samples",
             "pCR_mapped_samples", "DRFS_mapped_samples", "primary_TNBC_pCR_complete", "primary_TNBC_DRFS_complete"),
    n = c(
      ncol(expr_probe),
      nrow(clin),
      nrow(dat),
      sum(dat$primary_tnbc %in% TRUE, na.rm = TRUE),
      sum(dat$ER_HER2_negative_sensitivity %in% TRUE, na.rm = TRUE),
      sum(dat$basal_like_sensitivity %in% TRUE, na.rm = TRUE),
      sum(!is.na(dat$pCR_binary)),
      sum(!is.na(dat$drfs_event) & !is.na(dat$drfs_time_years)),
      sum(dat$primary_tnbc %in% TRUE & !is.na(dat$pCR_binary) & !is.na(dat$TNBC_BIR_favorable_score)),
      sum(dat$primary_tnbc %in% TRUE & !is.na(dat$drfs_event) & !is.na(dat$drfs_time_years) & !is.na(dat$TNBC_BIR_favorable_score))
    ),
    notes = c(
      "GSE25055 and GSE25065 processed series matrix expression columns.",
      "Parsed GEO sample phenotype rows.",
      "Rows retained after score-clinical join by GEO accession.",
      "Primary TNBC: ER-/PR-/HER2-.",
      "Sensitivity subset: ER-/HER2-, PR ignored.",
      "Sensitivity subset: PAM50 basal-like.",
      "Samples with mapped pCR/RD endpoint.",
      "Samples with mapped DRFS event and time.",
      "Complete primary TNBC pCR model cases.",
      "Complete primary TNBC DRFS Cox model cases."
    )
  )
  readr::write_tsv(count_flow, file.path(results_dir, "GSE25066_sample_count_flow.tsv"))

  if (count_flow$n[count_flow$step == "primary_TNBC_samples"] < 60 ||
      sum(dat$primary_tnbc %in% TRUE & dat$pCR_binary == 1, na.rm = TRUE) < 15) {
    message("GSE25066 primary TNBC subset is underpowered; models will be labeled but still run if estimable.")
  }

  pcr_res <- fit_pcr_models(dat)
  readr::write_tsv(pcr_res$model_tbl, file.path(results_dir, "GSE25066_pCR_models.tsv"))
  readr::write_tsv(pcr_res$auc_tbl, file.path(results_dir, "GSE25066_pCR_model_auc.tsv"))
  readr::write_tsv(pcr_res$lrt_tbl, file.path(results_dir, "GSE25066_pCR_model_lrt.tsv"))

  drfs_res <- fit_drfs_models(dat)
  readr::write_tsv(drfs_res$model_tbl, file.path(results_dir, "GSE25066_DRFS_models.tsv"))
  readr::write_tsv(drfs_res$ph_tbl, file.path(results_dir, "GSE25066_DRFS_PH_assumption.tsv"))

  make_figures(dat, pcr_res, drfs_res)
  message("GSE25066 external validation complete.")
}

main()
