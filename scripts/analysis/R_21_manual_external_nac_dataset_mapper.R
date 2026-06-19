#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(dplyr)
  library(org.Hs.eg.db)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(yaml)
})

source(file.path(getwd(), "scripts", "R", "utils_signature_scoring.R"))

set.seed(42)
root <- getwd()
cfg <- yaml::read_yaml(file.path(root, "config", "project_config.yaml"))
public_geo_root <- file.path(cfg$paths$public_data_root, "GEO")
results_dir <- file.path(root, cfg$paths$results)
processed_dir <- file.path(root, cfg$paths$processed_data)
reports_dir <- file.path(root, "reports")
logs_dir <- file.path(root, "logs")
figures_dir <- file.path(root, "figures")
dir.create(public_geo_root, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

priority_path <- file.path(root, "config", "external_validation_manual_dataset_priority.tsv")
priority <- readr::read_tsv(priority_path, show_col_types = FALSE)
sig_tbl <- readr::read_tsv(file.path(processed_dir, "TNBC_BIR_signatures_clean.tsv"), show_col_types = FALSE)
sig_tbl <- sig_tbl %>% mutate(gene = toupper(gene))

manual_accessions <- c("GSE41998", "GSE20194", "GSE22226", "GSE32646", "GSE163882")
core_model_scores <- c(
  "TNBC_BIR_favorable_score", "Immune_remodeling_score", "HLAII",
  "Proliferation", "APC_CXCL9_axis"
)

sha256_file <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NA_character_)
  out <- tryCatch(system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE), error = function(e) character())
  if (!length(out)) return(NA_character_)
  strsplit(out[[1]], "\\s+")[[1]][[1]]
}

one_line <- function(x) {
  x <- paste(as.character(x), collapse = " | ")
  x <- gsub("[\r\n\t]+", " ", x)
  gsub("\\s+", " ", trimws(x))
}

classify_download_error <- function(txt) {
  txt <- paste(txt, collapse = " ")
  dplyr::case_when(
    grepl("404|not found", txt, ignore.case = TRUE) ~ "http_404_or_missing",
    grepl("SSL|SYSCALL|handshake", txt, ignore.case = TRUE) ~ "ssl_or_dns_path",
    grepl("timed out|timeout|Could not resolve|Failed to connect", txt, ignore.case = TRUE) ~ "timeout_or_connectivity",
    TRUE ~ "download_failed_other"
  )
}

download_file <- function(url, dest, dataset, role) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && file.info(dest)$size > 0) {
    return(tibble(
      dataset = dataset, role = role, url = url, target_path = dest,
      status = "cached", error_class = NA_character_, message = "existing non-empty file",
      file_size = as.numeric(file.info(dest)$size), sha256 = sha256_file(dest),
      download_time = as.character(Sys.time())
    ))
  }
  args <- c(
    "--resolve", "ftp.ncbi.nlm.nih.gov:443:130.14.250.13",
    "-L", "--fail", "--retry", "3", "-C", "-", "-o", dest, url
  )
  out <- tryCatch(system2("curl", args, stdout = TRUE, stderr = TRUE), error = function(e) e)
  ok <- !inherits(out, "error") && file.exists(dest) && file.info(dest)$size > 0
  if (!ok && file.exists(dest) && file.info(dest)$size == 0) unlink(dest)
  msg <- if (inherits(out, "error")) conditionMessage(out) else one_line(out)
  tibble(
    dataset = dataset, role = role, url = url, target_path = dest,
    status = if (ok) "downloaded" else "failed",
    error_class = if (ok) NA_character_ else classify_download_error(msg),
    message = if (ok) "download completed" else msg,
    file_size = if (file.exists(dest)) as.numeric(file.info(dest)$size) else NA_real_,
    sha256 = sha256_file(dest),
    download_time = as.character(Sys.time())
  )
}

gse_bucket <- function(accession) {
  n <- as.integer(sub("^GSE", "", accession))
  paste0("GSE", floor(n / 1000), "nnn")
}

gpl_bucket <- function(platform) {
  n <- as.integer(sub("^GPL", "", platform))
  if (is.na(n)) stop("Unrecognized platform: ", platform)
  if (n < 1000) "GPLnnn" else paste0("GPL", floor(n / 1000), "nnn")
}

series_matrix_path <- function(accession) {
  file.path(public_geo_root, accession, paste0(accession, "_series_matrix.txt.gz"))
}

series_matrix_url <- function(accession) {
  paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/", gse_bucket(accession),
    "/", accession, "/matrix/", accession, "_series_matrix.txt.gz"
  )
}

download_series_matrix <- function(accession) {
  download_file(series_matrix_url(accession), series_matrix_path(accession), accession, "series_matrix")
}

platform_series_matrix_path <- function(accession, platform) {
  file.path(public_geo_root, accession, paste0(accession, "-", platform, "_series_matrix.txt.gz"))
}

platform_series_matrix_url <- function(accession, platform) {
  paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/", gse_bucket(accession),
    "/", accession, "/matrix/", accession, "-", platform, "_series_matrix.txt.gz"
  )
}

download_platform_series_matrix <- function(accession, platform) {
  download_file(
    platform_series_matrix_url(accession, platform),
    platform_series_matrix_path(accession, platform),
    accession,
    paste0("series_matrix_", platform)
  )
}

parse_matrix_line <- function(line) {
  read.delim(text = line, sep = "\t", header = FALSE, quote = "\"", comment.char = "",
             check.names = FALSE, stringsAsFactors = FALSE)
}

parse_pheno_from_series <- function(path) {
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
  if (is.null(geo)) stop("No !Sample_geo_accession in series matrix.")
  pheno <- tibble(sample_id = geo, title = if (!is.null(title)) title else geo)
  characteristic_rows <- parsed[vapply(parsed, function(x) identical(x[1, 1], "!Sample_characteristics_ch1"), logical(1))]
  for (row in characteristic_rows) {
    vals <- as.character(row[1, -1])
    keys <- ifelse(grepl(":", vals), sub(":.*$", "", vals), "characteristics_ch1")
    key <- names(sort(table(keys), decreasing = TRUE))[1]
    col <- key %>%
      stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
      stringr::str_replace_all("^_|_$", "") %>%
      tolower()
    value <- ifelse(grepl(":", vals), sub("^[^:]+:\\s*", "", vals), vals)
    if (col %in% names(pheno)) col <- make.unique(c(names(pheno), col))[length(names(pheno)) + 1]
    pheno[[col]] <- value
  }
  pheno
}

get_cached_pheno <- function(accession) {
  rds_path <- file.path(public_geo_root, accession, paste0(accession, "_GSEMatrix.rds"))
  if (file.exists(rds_path)) {
    obj <- readRDS(rds_path)
    if (inherits(obj, "light_gse")) return(obj$pheno)
    if (inherits(obj, "ExpressionSet")) return(Biobase::pData(obj) %>% as_tibble(rownames = "sample_id"))
    if (is.list(obj) && length(obj) > 0 && inherits(obj[[1]], "ExpressionSet")) {
      return(Biobase::pData(obj[[1]]) %>% as_tibble(rownames = "sample_id"))
    }
  }
  parse_pheno_from_series(series_matrix_path(accession))
}

parse_series_expression <- function(accession) {
  path <- series_matrix_path(accession)
  if (!file.exists(path)) stop("Missing series matrix: ", path)
  lines <- readLines(gzfile(path), warn = FALSE)
  begin <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end", lines)
  if (length(begin) != 1 || length(end) != 1 || end <= begin + 2) {
    stop("Series matrix expression table is absent or empty for ", accession)
  }
  tbl <- read.delim(
    text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"),
    sep = "\t", header = TRUE, quote = "\"", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE
  )
  id <- as.character(tbl[[1]])
  expr <- as.matrix(tbl[, -1, drop = FALSE])
  storage.mode(expr) <- "numeric"
  rownames(expr) <- id
  expr
}

parse_gse22226_platform <- function(platform) {
  accession <- "GSE22226"
  path <- platform_series_matrix_path(accession, platform)
  if (!file.exists(path)) stop("Missing GSE22226 platform matrix: ", path)
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
  pheno <- tibble(sample_id = geo, title = title, platform = platform)
  characteristic_rows <- parsed[vapply(parsed, function(x) identical(x[1, 1], "!Sample_characteristics_ch2"), logical(1))]
  for (j in seq_along(geo)) {
    vals <- as.character(vapply(characteristic_rows, function(row) {
      if (j + 1 <= ncol(row)) as.character(row[1, j + 1]) else ""
    }, character(1)))
    vals <- vals[nzchar(vals)]
    for (val in vals) {
      if (!grepl(":", val)) next
      key <- sub(":.*$", "", val)
      col <- key %>%
        stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
        stringr::str_replace_all("^_|_$", "") %>%
        tolower()
      value <- sub("^[^:]+:[[:space:]]*", "", val)
      pheno[j, col] <- value
    }
  }
  begin <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end", lines)
  if (length(begin) != 1 || length(end) != 1 || end <= begin + 2) {
    stop("GSE22226 platform expression table is absent or empty for ", platform)
  }
  tbl <- read.delim(
    text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"),
    sep = "\t", header = TRUE, quote = "\"", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE
  )
  expr <- as.matrix(tbl[, -1, drop = FALSE])
  storage.mode(expr) <- "numeric"
  rownames(expr) <- as.character(tbl[[1]])
  list(pheno = pheno, expr = expr)
}

read_gse22226_old_annotation <- function(platform) {
  accession <- "GSE22226"
  out_dir <- file.path(public_geo_root, accession, "annotation")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  raw_tar <- file.path(public_geo_root, accession, "GSE22226_RAW.tar")
  raw_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/", gse_bucket(accession),
    "/", accession, "/suppl/GSE22226_RAW.tar"
  )
  raw_dl <- download_file(raw_url, raw_tar, accession, "raw_tar_annotation_bundle")
  ann_file <- file.path(out_dir, paste0(platform, "_old_annotations.txt.gz"))
  if (!file.exists(ann_file) && file.exists(raw_tar)) {
    status <- system2("tar", c("-xf", raw_tar, "-C", out_dir), stdout = TRUE, stderr = TRUE)
    if (!file.exists(ann_file)) stop("Could not extract ", basename(ann_file), ": ", one_line(status))
  }
  lines <- readLines(gzfile(ann_file), warn = FALSE)
  start <- grep("^ID\t", lines)[1]
  if (is.na(start)) stop("No ID header in ", ann_file)
  tbl <- read.delim(
    text = paste(lines[start:length(lines)], collapse = "\n"),
    sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE,
    fill = TRUE, quote = ""
  )
  ann <- tibble(
    probe_id = as.character(tbl$ID),
    gene_symbol = vapply(tbl$GENE_SYMBOL, clean_symbol_one, character(1))
  ) %>%
    filter(!is.na(probe_id), probe_id != "", !is.na(gene_symbol), gene_symbol != "") %>%
    distinct(probe_id, .keep_all = TRUE)
  list(annotation = ann, download_qc = raw_dl, note = paste0(platform, " old annotation GENE_SYMBOL from GSE22226_RAW.tar"))
}

clean_symbol_one <- function(x) {
  x <- as.character(x)
  x <- gsub('"', "", x, fixed = TRUE)
  x <- gsub("\\s+", " ", x)
  parts <- unlist(strsplit(x, "///|//|;|,|\\|", perl = TRUE))
  parts <- trimws(parts)
  parts <- parts[parts != "" & !is.na(parts)]
  parts <- parts[!grepl("^(---|NA|N/A|null|unknown|LOC\\d+)$", parts, ignore.case = TRUE)]
  parts <- parts[grepl("^[A-Za-z][A-Za-z0-9_.-]*$", parts)]
  if (!length(parts)) return(NA_character_)
  toupper(parts[[1]])
}

read_platform_annotation <- function(accession, platform) {
  out_dir <- file.path(public_geo_root, accession)
  dest <- file.path(out_dir, paste0(platform, ".annot.gz"))
  url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/platforms/", gpl_bucket(platform), "/", platform, "/annot/", platform, ".annot.gz")
  dl <- download_file(url, dest, accession, paste0(platform, "_annotation"))
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    return(list(annotation = tibble(probe_id = character(), gene_symbol = character()), download_qc = dl,
                note = paste("annotation download failed for", platform)))
  }
  lines <- readLines(gzfile(dest), warn = FALSE)
  begin <- grep("^!platform_table_begin", lines)
  end <- grep("^!platform_table_end", lines)
  if (length(begin) == 1 && length(end) == 1 && end > begin + 1) {
    tbl <- read.delim(
      text = paste(lines[(begin + 1):(end - 1)], collapse = "\n"),
      sep = "\t", header = TRUE, quote = "\"", comment.char = "",
      check.names = FALSE, stringsAsFactors = FALSE
    )
  } else {
    tbl <- tryCatch(read.delim(gzfile(dest), sep = "\t", header = TRUE, quote = "\"",
                               comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE),
                    error = function(e) tibble())
  }
  if (!nrow(tbl)) {
    return(list(annotation = tibble(probe_id = character(), gene_symbol = character()), download_qc = dl,
                note = paste("annotation table empty for", platform)))
  }
  id_candidates <- c("ID", "ID_REF", "Probe Set ID", "Probe.Set.ID", "Reporter Identifier", "ProbeName", "Probe ID")
  id_col <- id_candidates[id_candidates %in% names(tbl)][1]
  if (is.na(id_col)) id_col <- names(tbl)[1]
  symbol_candidates <- names(tbl)[grepl("^gene symbol$|gene symbol|gene_symbol|symbol$", names(tbl), ignore.case = TRUE)]
  symbol_candidates <- symbol_candidates[!grepl("title|description|chromosome|unigene", symbol_candidates, ignore.case = TRUE)]
  symbol_col <- symbol_candidates[1]
  if (is.na(symbol_col)) {
    symbol_col <- names(tbl)[grepl("gene_assignment|gene assignment", names(tbl), ignore.case = TRUE)][1]
  }
  if (is.na(symbol_col)) {
    return(list(annotation = tibble(probe_id = as.character(tbl[[id_col]]), gene_symbol = NA_character_),
                download_qc = dl, note = paste("no gene symbol column detected for", platform)))
  }
  ann <- tibble(
    probe_id = as.character(tbl[[id_col]]),
    gene_symbol = vapply(tbl[[symbol_col]], clean_symbol_one, character(1))
  ) %>%
    filter(!is.na(probe_id), probe_id != "", !is.na(gene_symbol), gene_symbol != "") %>%
    distinct(probe_id, .keep_all = TRUE)
  list(annotation = ann, download_qc = dl, note = paste("symbol column:", symbol_col))
}

collapse_probe_expression <- function(expr, ann, rule = c("highest_iqr", "median")) {
  rule <- match.arg(rule)
  idx <- match(rownames(expr), ann$probe_id)
  symbol <- ann$gene_symbol[idx]
  keep <- !is.na(symbol) & symbol != ""
  expr <- expr[keep, , drop = FALSE]
  symbol <- symbol[keep]
  if (!nrow(expr)) stop("No probes mapped to gene symbols.")
  if (rule == "median") {
    df <- as.data.frame(expr) %>% mutate(gene_symbol = symbol)
    out <- df %>%
      group_by(gene_symbol) %>%
      summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)), .groups = "drop")
    mat <- as.matrix(out[, -1, drop = FALSE])
    rownames(mat) <- out$gene_symbol
    return(mat)
  }
  stat <- apply(expr, 1, IQR, na.rm = TRUE)
  keep_tbl <- tibble(gene_symbol = symbol, idx = seq_along(symbol), stat = stat) %>%
    group_by(gene_symbol) %>%
    slice_max(stat, n = 1, with_ties = FALSE) %>%
    ungroup()
  out <- expr[keep_tbl$idx, , drop = FALSE]
  rownames(out) <- keep_tbl$gene_symbol
  out
}

map_ensembl_to_symbol <- function(ids) {
  ens <- sub("\\..*$", "", ids)
  sym <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(ens),
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  out <- unname(sym[ens])
  toupper(out)
}

parse_gse163882_expression <- function(accession) {
  dest <- file.path(public_geo_root, accession, "GSE163882_all.data.tpms_222Samples.csv.gz")
  url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/", gse_bucket(accession),
    "/", accession, "/suppl/GSE163882_all.data.tpms_222Samples.csv.gz"
  )
  dl <- download_file(url, dest, accession, "supplementary_tpm_matrix")
  if (!file.exists(dest) || file.info(dest)$size == 0) stop("Missing GSE163882 supplementary TPM matrix.")
  tbl <- readr::read_csv(dest, show_col_types = FALSE, progress = FALSE)
  gene_id <- as.character(tbl[[1]])
  expr <- as.matrix(tbl[, -1, drop = FALSE])
  storage.mode(expr) <- "numeric"
  empty_cols <- colnames(expr) == "" | grepl("^\\.\\.\\.", colnames(expr))
  all_na_cols <- colSums(!is.na(expr)) == 0
  expr <- expr[, !(empty_cols | all_na_cols), drop = FALSE]
  symbol <- map_ensembl_to_symbol(gene_id)
  keep <- !is.na(symbol) & symbol != ""
  expr <- expr[keep, , drop = FALSE]
  symbol <- symbol[keep]
  stat <- apply(expr, 1, IQR, na.rm = TRUE)
  keep_tbl <- tibble(gene_symbol = symbol, idx = seq_along(symbol), stat = stat) %>%
    group_by(gene_symbol) %>%
    slice_max(stat, n = 1, with_ties = FALSE) %>%
    ungroup()
  out <- expr[keep_tbl$idx, , drop = FALSE]
  rownames(out) <- keep_tbl$gene_symbol
  list(expr = out, download_qc = dl, annotation_note = "org.Hs.eg.db ENSEMBL-to-SYMBOL mapping")
}

add_axis <- function(scores) {
  axis_parts <- c("HLAII", "CXCL9_FOLR2_macrophage", "DC_HLA_antigen_presentation")
  if (!all(axis_parts %in% names(scores))) {
    scores$APC_CXCL9_axis <- NA_real_
    return(scores)
  }
  zmat <- as.data.frame(lapply(scores[axis_parts], function(x) as.numeric(scale(x))))
  scores$APC_CXCL9_axis <- rowMeans(zmat, na.rm = TRUE)
  scores
}

norm_value <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[x %in% c("", "na", "nan", "unknown", "nd", "not available")] <- NA_character_
  x
}

is_neg <- function(x) norm_value(x) %in% c("negative", "neg", "n", "0", "no", "normal")
is_pos <- function(x) norm_value(x) %in% c("positive", "pos", "p", "1", "yes", "amplified")

map_pcr <- function(x) {
  v <- norm_value(x)
  out <- rep(NA_integer_, length(v))
  out[v %in% c("yes", "pcr", "pathologic complete response", "complete response")] <- 1L
  out[v %in% c("no", "rd", "ncr", "non-pcr", "non pcr", "residual disease")] <- 0L
  out
}

map_clinical <- function(accession, pheno) {
  pheno <- as_tibble(pheno)
  if (!"sample_id" %in% names(pheno)) pheno <- pheno %>% rownames_to_column("sample_id")
  get_col <- function(nm) if (nm %in% names(pheno)) pheno[[nm]] else rep(NA_character_, nrow(pheno))
  num_col <- function(nm) suppressWarnings(as.numeric(norm_value(get_col(nm))))
  if (accession == "GSE41998") {
    out <- pheno %>%
      transmute(
        dataset = accession, sample_id = sample_id, geo_sample_id = sample_id,
        title = title, tissue = NA_character_,
        endpoint_raw = as.character(pcr), pCR_binary = map_pcr(pcr),
        er_raw = as.character(er), pr_raw = as.character(pr), her2_raw = as.character(her2stat),
        er_negative = is_neg(er), pr_negative = is_neg(pr), her2_negative = is_neg(her2stat),
        cohort_note = paste0("treatment_arm=", treatment_arm)
      )
  } else if (accession == "GSE20194") {
    out <- pheno %>%
      transmute(
        dataset = accession, sample_id = sample_id, geo_sample_id = sample_id,
        title = title, tissue = tissue,
        endpoint_raw = as.character(pcr_vs_rd), pCR_binary = map_pcr(pcr_vs_rd),
        er_raw = as.character(er_status), pr_raw = as.character(pr_status), her2_raw = as.character(her2_status),
        er_negative = is_neg(er_status), pr_negative = is_neg(pr_status), her2_negative = is_neg(her2_status),
        cohort_note = paste0("tissue=", tissue)
      )
  } else if (accession == "GSE32646") {
    out <- pheno %>%
      transmute(
        dataset = accession, sample_id = sample_id, geo_sample_id = sample_id,
        title = title, tissue = tissue,
        endpoint_raw = as.character(pathologic_response_pcr_ncr), pCR_binary = map_pcr(pathologic_response_pcr_ncr),
        er_raw = as.character(er_status_ihc), pr_raw = as.character(pr_status_ihc), her2_raw = as.character(her2_status_fish),
        er_negative = is_neg(er_status_ihc), pr_negative = is_neg(pr_status_ihc), her2_negative = is_neg(her2_status_fish),
        cohort_note = paste0("stage=", clinical_stage)
      )
  } else if (accession == "GSE163882") {
    sample_key <- sub(":.*$", "", as.character(pheno$title))
    out <- pheno %>%
      mutate(sample_key = sample_key) %>%
      transmute(
        dataset = accession, sample_id = sample_key, geo_sample_id = sample_id,
        title = title, tissue = tissue,
        endpoint_raw = as.character(response_to_nac), pCR_binary = map_pcr(response_to_nac),
        er_raw = as.character(estrogen_receptor_status), pr_raw = as.character(progesterone_receptor_status),
        her2_raw = as.character(her2_receptor_status),
        er_negative = is_neg(estrogen_receptor_status), pr_negative = is_neg(progesterone_receptor_status),
        her2_negative = is_neg(her2_receptor_status),
        cohort_note = paste0("stage=", breast_cancer_stage)
      )
  } else if (accession == "GSE22226") {
    rfs_time_raw <- get_col("relapse_free_survival_time_time_from_chemo_start_date_until_earliest")
    rfs_time <- suppressWarnings(as.numeric(sub(".*days\\):[[:space:]]*", "", as.character(rfs_time_raw))))
    out <- pheno %>%
      transmute(
        dataset = accession, sample_id = sample_id, geo_sample_id = sample_id,
        title = title, tissue = "pretreatment breast tumor biopsy",
        platform = platform,
        endpoint_raw = as.character(get_col("pathological_complete_response_pcr")),
        pCR_binary = map_pcr(get_col("pathological_complete_response_pcr")),
        er_raw = as.character(get_col("er_0_negative_1_positive")),
        pr_raw = as.character(get_col("pgr_0_negative_1_positive")),
        her2_raw = as.character(get_col("her2_0_negative_1_positive")),
        er_negative = num_col("er_0_negative_1_positive") == 0,
        pr_negative = num_col("pgr_0_negative_1_positive") == 0,
        her2_negative = num_col("her2_0_negative_1_positive") == 0,
        pam50_basal_like = norm_value(get_col("intrinsic_subtype_by_pam50")) == "basal-like",
        rfs_time_days = rfs_time,
        rfs_event = suppressWarnings(as.integer(norm_value(get_col("relapse_free_survival_indicator_1_event_local_or_distant_progression_or_death_0_censor_at_last_follow_up")))),
        cohort_note = paste0("platform=", platform, "; pam50=", get_col("intrinsic_subtype_by_pam50"))
      )
  } else {
    stop("No clinical mapper for ", accession)
  }
  out %>%
    mutate(
      pCR_mapped = !is.na(pCR_binary),
      primary_tnbc = er_negative & pr_negative & her2_negative,
      er_her2_negative = er_negative & her2_negative,
      pam50_basal_like = if ("pam50_basal_like" %in% names(.)) pam50_basal_like else NA,
      subtype_mapping_rule = "clinical receptor status from GEO sample characteristics"
    )
}

score_dataset <- function(accession, expr_gene, clinical, expression_source, annotation_note) {
  scored <- score_signatures(expr_gene, sig_tbl, min_coverage = 0.70, duplicate_rule = "highest_mad")
  scores <- add_axis(scored$scores)
  score_cols <- intersect(c(names(scored$scores), "APC_CXCL9_axis"), names(scores))
  joined <- clinical %>%
    inner_join(scores %>% select(sample_id, all_of(setdiff(score_cols, "sample_id"))), by = "sample_id") %>%
    mutate(expression_source = expression_source, annotation_note = annotation_note)
  coverage <- scored$coverage %>%
    mutate(
      dataset = accession,
      expression_source = expression_source,
      annotation_note = annotation_note,
      .before = 1
    )
  axis_cov <- coverage %>%
    filter(signature %in% c("HLAII", "CXCL9_FOLR2_macrophage", "DC_HLA_antigen_presentation")) %>%
    summarise(
      dataset = accession,
      expression_source = .env$expression_source,
      annotation_note = .env$annotation_note,
      signature = "APC_CXCL9_axis",
      n_signature_genes = sum(n_signature_genes),
      n_genes_found = sum(n_genes_found),
      coverage = min(coverage, na.rm = TRUE),
      missing_genes = paste(missing_genes, collapse = ";")
    )
  coverage <- bind_rows(coverage, axis_cov)
  list(scores_clinical = joined, coverage = coverage)
}

make_qc_row <- function(accession, clinical, joined, coverage, expression_features, expression_samples,
                        expression_source, note) {
  primary <- joined %>% filter(primary_tnbc, pCR_mapped)
  erher2 <- joined %>% filter(er_her2_negative, pCR_mapped)
  modeled_cov <- coverage %>%
    filter(signature %in% c("Immune_remodeling_score", "HLAII", "Proliferation", "APC_CXCL9_axis")) %>%
    summarise(min_cov = min(coverage, na.rm = TRUE)) %>%
    pull(min_cov)
  if (!is.finite(modeled_cov)) modeled_cov <- NA_real_
  n_primary <- nrow(primary)
  p_primary <- sum(primary$pCR_binary == 1, na.rm = TRUE)
  rd_primary <- sum(primary$pCR_binary == 0, na.rm = TRUE)
  decision <- case_when(
    is.na(modeled_cov) | modeled_cov < 0.70 ~ "blocker_low_signature_coverage",
    n_primary >= 25 & p_primary >= 8 & rd_primary >= 8 ~ "eligible_for_external_validation_modeling",
    TRUE ~ "not_eligible_underpowered_or_unmapped"
  )
  reason <- case_when(
    decision == "eligible_for_external_validation_modeling" ~ "Primary TNBC pCR subset passes n/event and signature coverage criteria.",
    decision == "blocker_low_signature_coverage" ~ "Modeled score signature coverage is below 70%.",
    TRUE ~ "Primary TNBC pCR subset does not meet n>=25 and pCR/RD>=8 eligibility criteria."
  )
  tibble(
    dataset = accession,
    priority = priority$priority[match(accession, priority$accession)],
    expression_source = expression_source,
    expression_features = expression_features,
    expression_samples = expression_samples,
    clinical_samples = nrow(clinical),
    joined_samples = nrow(joined),
    primary_tnbc_n = n_primary,
    primary_tnbc_pcr = p_primary,
    primary_tnbc_rd = rd_primary,
    er_her2_negative_n = nrow(erher2),
    er_her2_negative_pcr = sum(erher2$pCR_binary == 1, na.rm = TRUE),
    er_her2_negative_rd = sum(erher2$pCR_binary == 0, na.rm = TRUE),
    endpoint_mapped_n = sum(joined$pCR_mapped, na.rm = TRUE),
    subtype_mapped_n = sum(!is.na(joined$primary_tnbc), na.rm = TRUE),
    min_modeled_signature_coverage = modeled_cov,
    coverage_pass = isTRUE(modeled_cov >= 0.70),
    decision = decision,
    reason = reason,
    notes = note
  )
}

process_affy_dataset <- function(accession, platform) {
  dl_series <- download_series_matrix(accession)
  pheno <- get_cached_pheno(accession)
  expr <- parse_series_expression(accession)
  ann_obj <- read_platform_annotation(accession, platform)
  expr_gene <- collapse_probe_expression(expr, ann_obj$annotation, rule = "highest_iqr")
  clinical <- map_clinical(accession, pheno)
  scored <- score_dataset(accession, expr_gene, clinical, "GEO_series_matrix_highest_IQR_probe_collapse", ann_obj$note)
  qc <- make_qc_row(
    accession, clinical, scored$scores_clinical, scored$coverage,
    nrow(expr_gene), ncol(expr_gene),
    "GEO_series_matrix_highest_IQR_probe_collapse", ann_obj$note
  )
  list(
    scores_clinical = scored$scores_clinical,
    coverage = scored$coverage,
    qc = qc,
    download_qc = bind_rows(dl_series, ann_obj$download_qc)
  )
}

process_gse163882 <- function() {
  accession <- "GSE163882"
  dl_series <- download_series_matrix(accession)
  pheno <- get_cached_pheno(accession)
  parsed <- parse_gse163882_expression(accession)
  clinical <- map_clinical(accession, pheno)
  scored <- score_dataset(accession, parsed$expr, clinical, "supplementary_TPM_ENSEMBL_orgHs_symbol", parsed$annotation_note)
  qc <- make_qc_row(
    accession, clinical, scored$scores_clinical, scored$coverage,
    nrow(parsed$expr), ncol(parsed$expr),
    "supplementary_TPM_ENSEMBL_orgHs_symbol", parsed$annotation_note
  )
  list(
    scores_clinical = scored$scores_clinical,
    coverage = scored$coverage,
    qc = qc,
    download_qc = bind_rows(dl_series, parsed$download_qc)
  )
}

process_gse22226 <- function() {
  accession <- "GSE22226"
  out_dir <- file.path(public_geo_root, accession)
  filelist_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/", gse_bucket(accession),
    "/", accession, "/suppl/filelist.txt"
  )
  fl_dl <- download_file(filelist_url, file.path(out_dir, "filelist.txt"), accession, "supplementary_filelist")
  platform_results <- lapply(c("GPL1708", "GPL4133"), function(platform) {
    dl_series <- download_platform_series_matrix(accession, platform)
    parsed <- parse_gse22226_platform(platform)
    ann_obj <- read_gse22226_old_annotation(platform)
    expr_gene <- collapse_probe_expression(parsed$expr, ann_obj$annotation, rule = "highest_iqr")
    clinical <- map_clinical(accession, parsed$pheno)
    source_label <- paste0("GEO_platform_series_matrix_", platform, "_highest_IQR_probe_collapse")
    scored <- score_dataset(accession, expr_gene, clinical, source_label, ann_obj$note)
    list(
      scores_clinical = scored$scores_clinical,
      coverage = scored$coverage %>% mutate(platform = platform, .after = dataset),
      clinical = clinical,
      expression_features = nrow(expr_gene),
      expression_samples = ncol(expr_gene),
      download_qc = bind_rows(dl_series, ann_obj$download_qc)
    )
  })
  scores_clinical <- bind_rows(lapply(platform_results, `[[`, "scores_clinical"))
  coverage <- bind_rows(lapply(platform_results, `[[`, "coverage"))
  clinical <- bind_rows(lapply(platform_results, `[[`, "clinical"))
  qc <- make_qc_row(
    accession, clinical, scores_clinical, coverage,
    sum(vapply(platform_results, `[[`, numeric(1), "expression_features")),
    sum(vapply(platform_results, `[[`, numeric(1), "expression_samples")),
    "GEO_platform_series_matrix_GPL1708_GPL4133_highest_IQR_probe_collapse",
    "GSE22226 platform-specific series matrices plus old GPL annotation bundle; model should adjust for platform."
  ) %>%
    mutate(
      notes = paste(notes, "pCR events are exactly at the Round 4 minimum threshold; use platform-adjusted model and mark underpowered.", sep = "; ")
    )
  if (file.exists(file.path(logs_dir, "GSE22226_external_validation_manual_mapping_blocker.md"))) {
    unlink(file.path(logs_dir, "GSE22226_external_validation_manual_mapping_blocker.md"))
  }
  primary <- scores_clinical %>% filter(primary_tnbc, pCR_mapped)
  primary_n <- nrow(primary)
  primary_pcr <- sum(primary$pCR_binary == 1, na.rm = TRUE)
  primary_rd <- sum(primary$pCR_binary == 0, na.rm = TRUE)
  writeLines(c(
    "# GSE22226 Round 4 Mapping Recovery",
    "",
    paste("Generated:", Sys.time()),
    "",
    "GSE22226 was recovered using platform-specific GEO series matrices:",
    "- GSE22226-GPL1708_series_matrix.txt.gz",
    "- GSE22226-GPL4133_series_matrix.txt.gz",
    "",
    "Clinical metadata were parsed from per-sample `!Sample_characteristics_ch2` key-value entries, because the dual-channel reference design stores tumor/sample metadata on channel 2.",
    "",
    paste0(
      "The primary TNBC subset passes the minimum Round 4 threshold but is underpowered: n=",
      primary_n, ", pCR=", primary_pcr, ", RD=", primary_rd, "."
    ),
    "",
    "Modeling uses a platform fixed effect where possible."
  ), file.path(logs_dir, "GSE22226_external_validation_mapping_recovery.md"))
  list(
    scores_clinical = scores_clinical,
    coverage = coverage,
    qc = qc,
    download_qc = bind_rows(fl_dl, bind_rows(lapply(platform_results, `[[`, "download_qc")))
  )
}

low_priority_qc_rows <- function() {
  priority %>%
    filter(!accession %in% c("GSE25066", manual_accessions)) %>%
    transmute(
      dataset = accession,
      priority = priority,
      expression_source = NA_character_,
      expression_features = NA_integer_,
      expression_samples = NA_integer_,
      clinical_samples = NA_integer_,
      joined_samples = NA_integer_,
      primary_tnbc_n = NA_integer_,
      primary_tnbc_pcr = NA_integer_,
      primary_tnbc_rd = NA_integer_,
      er_her2_negative_n = NA_integer_,
      er_her2_negative_pcr = NA_integer_,
      er_her2_negative_rd = NA_integer_,
      endpoint_mapped_n = NA_integer_,
      subtype_mapped_n = NA_integer_,
      min_modeled_signature_coverage = NA_real_,
      coverage_pass = FALSE,
      decision = "deferred_lower_priority_external_validation",
      reason = "Round 4 manual mapping prioritized GSE41998, GSE20194, GSE22226, GSE32646, and GSE163882 as instructed; lower-priority datasets are documented but not modeled unless metadata are clearly usable.",
      notes = why_include
    )
}

gse25066_qc_row <- function() {
  tibble(
    dataset = "GSE25066",
    priority = priority$priority[match("GSE25066", priority$accession)],
    expression_source = "handled_by_scripts_19_and_20",
    expression_features = NA_integer_,
    expression_samples = NA_integer_,
    clinical_samples = NA_integer_,
    joined_samples = NA_integer_,
    primary_tnbc_n = NA_integer_,
    primary_tnbc_pcr = NA_integer_,
    primary_tnbc_rd = NA_integer_,
    er_her2_negative_n = NA_integer_,
    er_her2_negative_pcr = NA_integer_,
    er_her2_negative_rd = NA_integer_,
    endpoint_mapped_n = NA_integer_,
    subtype_mapped_n = NA_integer_,
    min_modeled_signature_coverage = NA_real_,
    coverage_pass = NA,
    decision = "handled_as_decisive_negative_technical_audit",
    reason = "GSE25066 processed and RMA analyses are generated by scripts/R/19 and scripts/R/20 and are reported separately.",
    notes = "Not re-entered into manual rescue mapping."
  )
}

main <- function() {
  results <- list()
  results[["GSE41998"]] <- tryCatch(process_affy_dataset("GSE41998", "GPL571"), error = function(e) e)
  results[["GSE20194"]] <- tryCatch(process_affy_dataset("GSE20194", "GPL96"), error = function(e) e)
  results[["GSE22226"]] <- tryCatch(process_gse22226(), error = function(e) e)
  results[["GSE32646"]] <- tryCatch(process_affy_dataset("GSE32646", "GPL570"), error = function(e) e)
  results[["GSE163882"]] <- tryCatch(process_gse163882(), error = function(e) e)

  qc_rows <- list(gse25066_qc_row())
  score_rows <- list()
  coverage_rows <- list()
  download_rows <- list()

  for (accession in names(results)) {
    obj <- results[[accession]]
    if (inherits(obj, "error")) {
      msg <- one_line(conditionMessage(obj))
      writeLines(c(
        paste0("# ", accession, " Round 4 Manual Mapping Blocker"),
        "",
        paste("Generated:", Sys.time()),
        "",
        paste("Error:", msg),
        "",
        "No model was run for this dataset. Re-run scripts/R/21 after resolving the blocker."
      ), file.path(logs_dir, paste0(accession, "_external_validation_manual_mapping_blocker.md")))
      qc_rows[[length(qc_rows) + 1]] <- tibble(
        dataset = accession,
        priority = priority$priority[match(accession, priority$accession)],
        expression_source = NA_character_,
        expression_features = NA_integer_,
        expression_samples = NA_integer_,
        clinical_samples = NA_integer_,
        joined_samples = NA_integer_,
        primary_tnbc_n = NA_integer_,
        primary_tnbc_pcr = NA_integer_,
        primary_tnbc_rd = NA_integer_,
        er_her2_negative_n = NA_integer_,
        er_her2_negative_pcr = NA_integer_,
        er_her2_negative_rd = NA_integer_,
        endpoint_mapped_n = NA_integer_,
        subtype_mapped_n = NA_integer_,
        min_modeled_signature_coverage = NA_real_,
        coverage_pass = FALSE,
        decision = "blocker_mapping_error",
        reason = msg,
        notes = paste0("See logs/", accession, "_external_validation_manual_mapping_blocker.md")
      )
    } else {
      qc_rows[[length(qc_rows) + 1]] <- obj$qc
      if (nrow(obj$scores_clinical)) score_rows[[length(score_rows) + 1]] <- obj$scores_clinical
      if (nrow(obj$coverage)) coverage_rows[[length(coverage_rows) + 1]] <- obj$coverage
      if (nrow(obj$download_qc)) download_rows[[length(download_rows) + 1]] <- obj$download_qc
    }
  }
  qc <- bind_rows(qc_rows, low_priority_qc_rows()) %>% arrange(priority, dataset)
  scores_clinical <- bind_rows(score_rows)
  coverage <- bind_rows(coverage_rows)
  download_qc <- bind_rows(download_rows)

  readr::write_tsv(qc, file.path(results_dir, "external_validation_external_dataset_mapping_qc.tsv"))
  readr::write_tsv(scores_clinical, file.path(processed_dir, "external_validation_external_NAC_scores_clinical.tsv"))
  readr::write_tsv(coverage, file.path(results_dir, "external_validation_external_signature_coverage.tsv"))
  readr::write_tsv(download_qc, file.path(results_dir, "external_validation_external_download_qc.tsv"))

  eligibility_summary <- qc %>%
    select(dataset, priority, decision, primary_tnbc_n, primary_tnbc_pcr, primary_tnbc_rd,
           min_modeled_signature_coverage, reason)
  report_lines <- c(
    "# External NAC Manual Mapping Report",
    "",
    paste("Generated:", Sys.time()),
    "",
    "Round 4 manual mapping was run before any external NAC modeling. The script used only locked TNBC-BIR signatures, fixed receptor/endpoint rules, and platform/gene annotation collapse; no outcome-derived cutoff or gene rescue was used.",
    "",
    "## Eligibility Summary",
    "",
    paste(capture.output(print(eligibility_summary, n = Inf)), collapse = "\n"),
    "",
    "## Notes",
    "",
    "- GSE25066 remains the decisive negative cohort and is handled by scripts 19/20, not by this rescue mapper.",
    "- GSE22226 is recovered via platform-specific series matrices when available; if recovery fails, the mapping QC table records the blocker.",
    "- Eligible rows are passed unchanged to scripts/R/22_manual_external_nac_validation_and_meta.R."
  )
  writeLines(report_lines, file.path(reports_dir, "external_NAC_manual_mapping_report.md"))

  message("Round 4 external NAC manual mapping complete.")
  invisible(TRUE)
}

main()
