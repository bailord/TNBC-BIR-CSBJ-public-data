#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(yaml)
})

source("scripts/R/utils_signature_scoring.R")

root <- getwd()
cfg <- yaml::read_yaml(file.path(root, "config", "project_config.yaml"))
sig_tbl <- readr::read_tsv(file.path(root, "data/processed/TNBC_BIR_signatures_clean.tsv"), show_col_types = FALSE)
if (nrow(sig_tbl) == 0) stop("Run scripts/R/01_prepare_signatures.R first.")

normalize_ispy2_id <- function(x) {
  x <- as.character(x)
  x <- sub("^ISPY2[_-]?", "", x, ignore.case = TRUE)
  x <- sub("-GPL[0-9]+$", "", x, ignore.case = TRUE)
  x <- sub("^X(?=\\d)", "", x, perl = TRUE)
  x <- trimws(x)
  x
}

parse_matrix_line <- function(line) {
  read.delim(
    text = line,
    sep = "\t",
    header = FALSE,
    quote = "\"",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

parse_geo_series_matrix <- function(path) {
  message("Parsing GEO series matrix: ", path)
  lines <- readLines(gzfile(path), warn = FALSE)
  sample_lines <- lines[grepl("^!Sample_", lines)]
  parsed <- lapply(sample_lines, parse_matrix_line)

  get_values <- function(key) {
    hits <- parsed[vapply(parsed, function(x) identical(x[1, 1], key), logical(1))]
    if (length(hits) == 0) return(NULL)
    as.character(hits[[1]][1, -1])
  }

  title <- get_values("!Sample_title")
  geo <- get_values("!Sample_geo_accession")
  if (is.null(title) || is.null(geo)) {
    stop("Could not find !Sample_title and !Sample_geo_accession in series matrix.")
  }

  pheno <- tibble(
    geo_accession = geo,
    sample_title = title,
    sample_id = normalize_ispy2_id(title)
  )

  characteristic_rows <- parsed[vapply(parsed, function(x) identical(x[1, 1], "!Sample_characteristics_ch1"), logical(1))]
  for (row in characteristic_rows) {
    vals <- as.character(row[1, -1])
    keys <- sub(":.*$", "", vals)
    key <- names(sort(table(keys), decreasing = TRUE))[1]
    col <- key %>%
      str_replace_all("[^A-Za-z0-9]+", "_") %>%
      str_replace_all("^_|_$", "")
    value <- ifelse(grepl(":", vals), sub("^[^:]+:\\s*", "", vals), vals)
    if (col %in% names(pheno)) col <- make.unique(c(names(pheno), col))[length(names(pheno)) + 1]
    pheno[[col]] <- value
  }

  if ("patient_id" %in% names(pheno)) {
    pheno$sample_id <- normalize_ispy2_id(pheno$patient_id)
  }
  pheno
}

read_ispy2_expression <- function(path) {
  message("Reading GSE194040 gene-level expression: ", path)
  expr_tbl <- utils::read.delim(
    gzfile(path),
    sep = "\t",
    header = FALSE,
    fill = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (nrow(expr_tbl) < 2 || ncol(expr_tbl) < 100) stop("Expression file has too few rows/columns: ", path)

  header <- as.character(unlist(expr_tbl[1, ], use.names = FALSE))
  body <- expr_tbl[-1, , drop = FALSE]
  gene <- toupper(trimws(as.character(body[[1]])))
  expr_body <- body[, -1, drop = FALSE]
  header <- header[seq_len(ncol(expr_body))]
  valid_samples <- !is.na(header) & header != ""
  expr_body <- expr_body[, valid_samples, drop = FALSE]
  header <- header[valid_samples]

  raw_sample_names <- normalize_ispy2_id(header)
  mat <- as.matrix(as.data.frame(lapply(expr_body, function(x) suppressWarnings(as.numeric(x))), check.names = FALSE))
  rownames(mat) <- gene
  colnames(mat) <- raw_sample_names
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  raw_n <- ncol(mat)
  dup_ids <- unique(colnames(mat)[duplicated(colnames(mat))])
  if (length(dup_ids) > 0) {
    message("Collapsing duplicate expression sample IDs by per-gene mean: ", paste(dup_ids, collapse = ", "))
    collapsed <- lapply(unique(colnames(mat)), function(sample_id) {
      idx <- which(colnames(mat) == sample_id)
      if (length(idx) == 1) mat[, idx] else rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
    })
    mat <- do.call(cbind, collapsed)
    colnames(mat) <- unique(raw_sample_names)
    rownames(mat) <- gene
  }
  attr(mat, "raw_sample_n") <- raw_n
  attr(mat, "duplicate_sample_ids") <- dup_ids
  mat
}

write_clean_expression <- function(expr, out_path) {
  out <- as.data.frame(expr, check.names = FALSE) %>%
    rownames_to_column("gene")
  readr::write_tsv(out, out_path)
}

score_gse194040 <- function() {
  gse <- "GSE194040"
  raw_dir <- file.path(root, cfg$paths$raw_data, gse)
  expr_path <- file.path(raw_dir, "GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt.gz")
  series_paths <- list.files(raw_dir, pattern = "GSE194040.*series_matrix[.]txt[.]gz$", full.names = TRUE)
  if (!file.exists(expr_path)) stop("Missing expression file. Run scripts/R/02_download_geo_bulk.R first: ", expr_path)
  if (length(series_paths) == 0) stop("Missing series matrix files. Run scripts/R/02_download_geo_bulk.R first in: ", raw_dir)

  expr <- read_ispy2_expression(expr_path)
  expr_raw_n <- attr(expr, "raw_sample_n")
  expr_dup_ids <- attr(expr, "duplicate_sample_ids")

  pheno_all <- bind_rows(lapply(series_paths, function(path) {
    parse_geo_series_matrix(path) %>% mutate(series_matrix_file = basename(path))
  }))
  pheno_dup_ids <- unique(pheno_all$sample_id[duplicated(pheno_all$sample_id)])
  pheno <- pheno_all %>%
    distinct(sample_id, .keep_all = TRUE)

  common <- intersect(colnames(expr), pheno$sample_id)
  if (length(common) < 100) {
    stop("Expression/clinical overlap is too small: ", length(common), " samples.")
  }
  expr <- expr[, common, drop = FALSE]
  pheno <- pheno %>% filter(sample_id %in% common) %>% arrange(match(sample_id, common))

  res <- score_signatures(
    expr,
    sig_tbl = sig_tbl,
    min_coverage = cfg$score_settings$min_signature_coverage,
    duplicate_rule = cfg$score_settings$duplicate_gene_rule
  )
  scores <- assign_bir_quadrants(res$scores)
  coverage <- res$coverage %>% mutate(dataset = gse, .before = 1)

  processed_dir <- file.path(root, cfg$paths$processed_data)
  results_dir <- file.path(root, cfg$paths$results)
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  write_clean_expression(expr, file.path(processed_dir, "GSE194040_expression.tsv.gz"))
  readr::write_tsv(pheno, file.path(processed_dir, "GSE194040_pheno_raw.tsv"))
  readr::write_tsv(scores, file.path(results_dir, "GSE194040_TNBC_BIR_scores.tsv"))
  readr::write_tsv(coverage, file.path(results_dir, "GSE194040_signature_coverage.tsv"))

  duplicate_qc <- tibble(
    dataset = gse,
    sample_id = union(expr_dup_ids, pheno_dup_ids),
    duplicated_in_expression = sample_id %in% expr_dup_ids,
    duplicated_in_clinical_series_matrix = sample_id %in% pheno_dup_ids,
    resolution = "Collapsed to one patient-level sample; duplicate expression columns were averaged by gene and duplicate clinical rows were de-duplicated by sample_id."
  )
  readr::write_tsv(duplicate_qc, file.path(results_dir, "GSE194040_duplicate_sample_id_qc.tsv"))

  count_flow <- tibble(
    dataset = gse,
    step = c(
      "expression_raw_columns",
      "expression_unique_samples_after_duplicate_collapse",
      "clinical_series_matrix_raw_rows",
      "clinical_unique_samples_after_duplicate_collapse",
      "expression_clinical_overlap"
    ),
    n = c(expr_raw_n, ncol(expr), nrow(pheno_all), nrow(pheno), length(common)),
    notes = c(
      "Raw sample columns in gene-level supplementary expression file after ID normalization.",
      "Unique patient-level expression samples after duplicate sample IDs were collapsed.",
      "Raw GEO series-matrix rows after parsing both platform matrix files.",
      "Unique patient-level clinical samples after duplicate sample IDs were collapsed.",
      "Samples retained for scoring after expression-clinical intersection."
    )
  )
  readr::write_tsv(count_flow, file.path(results_dir, "GSE194040_sample_count_flow_initial.tsv"))

  saveRDS(
    list(expr = expr, pheno = pheno, scores = scores, coverage = coverage, count_flow = count_flow),
    file.path(processed_dir, "GSE194040_scored_bulk.rds")
  )
  message("Scored ", gse, ": ", ncol(expr), " matched expression/clinical samples.")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- "GSE194040"

for (gse in args) {
  if (gse == "GSE194040") {
    score_gse194040()
  } else {
    warning("No direct scorer implemented yet for ", gse, ". Skipping.")
  }
}
