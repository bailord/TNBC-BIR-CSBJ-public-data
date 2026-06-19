#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(readr)
  library(tibble)
  library(tidyr)
  library(yaml)
})

root <- getwd()
cfg <- yaml::read_yaml(file.path(root, "config", "project_config.yaml"))
public_geo_root <- file.path(cfg$paths$public_data_root, "GEO")
results_dir <- file.path(root, cfg$paths$results)
processed_dir <- file.path(root, cfg$paths$processed_data)
figures_dir <- file.path(root, "figures")
reports_dir <- file.path(root, "reports")
logs_dir <- file.path(root, "logs")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

pkg_status <- function(stage) {
  pkgs <- c("SingleCellExperiment", "SummarizedExperiment", "Matrix", "SeuratObject", "Seurat")
  tibble(stage = stage, package = pkgs, available = vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)),
         check_time = as.character(Sys.time()))
}

ensure_single_cell_experiment <- function() {
  if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    return(tibble(package = "SingleCellExperiment", attempted = FALSE, success = TRUE,
                  message = "already installed", log_path = NA_character_))
  }
  log_path <- file.path(logs_dir, "external_validation_install_SingleCellExperiment.log")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    writeLines("BiocManager is unavailable; cannot install SingleCellExperiment.", log_path)
    return(tibble(package = "SingleCellExperiment", attempted = FALSE, success = FALSE,
                  message = "BiocManager unavailable", log_path = log_path))
  }
  out <- tryCatch(
    utils::capture.output(BiocManager::install("SingleCellExperiment", ask = FALSE, update = FALSE), type = "message"),
    error = function(e) e
  )
  if (inherits(out, "error")) writeLines(conditionMessage(out), log_path) else writeLines(out, log_path)
  tibble(package = "SingleCellExperiment", attempted = TRUE,
         success = requireNamespace("SingleCellExperiment", quietly = TRUE),
         message = if (requireNamespace("SingleCellExperiment", quietly = TRUE)) "install succeeded" else "install failed",
         log_path = log_path)
}

one_line <- function(x) {
  x <- paste(as.character(x), collapse = " | ")
  x <- gsub("[\r\n\t]+", " ", x)
  gsub("\\s+", " ", trimws(x))
}

decompress_double_gzip <- function(src, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && file.info(dest)$size > 0) return("cached")
  cmd <- paste("gzip -cd", shQuote(src), "| gzip -cd >", shQuote(dest))
  status <- system2("bash", c("-lc", cmd), stdout = TRUE, stderr = TRUE)
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    stop("double-gzip decompression failed: ", one_line(status))
  }
  "decompressed"
}

read_double_gzip_rds <- function(src) {
  con1 <- gzfile(src, "rb")
  con2 <- gzcon(con1)
  on.exit(close(con2), add = TRUE)
  readRDS(con2)
}

sig_tbl <- readr::read_tsv(file.path(processed_dir, "TNBC_BIR_signatures_clean.tsv"), show_col_types = FALSE) %>%
  mutate(gene = toupper(gene))
module_names <- c(
  "HLAII", "CXCL9_FOLR2_macrophage", "DC_HLA_antigen_presentation",
  "Mast_cell", "Tfh_TLS_like", "CD8_Tsem_memory",
  "Cytolytic", "Proliferation", "Suppressive_myeloid", "Basoluminal_plasticity"
)
sig_list <- split(toupper(sig_tbl$gene), sig_tbl$signature)
sig_list <- sig_list[names(sig_list) %in% module_names]

score_modules_by_cell <- function(counts, genes_upper, sig_list) {
  lib <- Matrix::colSums(counts)
  lib[lib <= 0 | is.na(lib)] <- 1
  out <- list()
  cov <- list()
  for (sig in names(sig_list)) {
    genes <- unique(sig_list[[sig]])
    present <- intersect(genes, genes_upper)
    idx <- match(present, genes_upper)
    cov[[sig]] <- tibble(
      signature = sig,
      n_signature_genes = length(genes),
      n_genes_found = length(present),
      coverage = length(present) / length(genes),
      missing_genes = paste(setdiff(genes, genes_upper), collapse = ";")
    )
    if (!length(idx)) {
      out[[sig]] <- rep(NA_real_, ncol(counts))
      next
    }
    sub <- counts[idx, , drop = FALSE]
    norm <- t(t(sub) / lib * 10000)
    out[[sig]] <- Matrix::colMeans(log1p(norm), na.rm = TRUE)
  }
  list(scores = as_tibble(out), coverage = bind_rows(cov))
}

safe_z <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (all(is.na(x)) || !is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

add_composites <- function(df) {
  axis_parts <- c("HLAII", "CXCL9_FOLR2_macrophage", "DC_HLA_antigen_presentation")
  immune_parts <- c("Mast_cell", "Tfh_TLS_like", "CD8_Tsem_memory", "CXCL9_FOLR2_macrophage", "DC_HLA_antigen_presentation")
  if (all(axis_parts %in% names(df))) {
    z <- as.data.frame(lapply(df[axis_parts], safe_z))
    df$APC_CXCL9_axis <- rowMeans(z, na.rm = TRUE)
  } else {
    df$APC_CXCL9_axis <- NA_real_
  }
  if (all(immune_parts %in% names(df))) {
    z <- as.data.frame(lapply(df[immune_parts], safe_z))
    df$Immune_remodeling_score <- rowMeans(z, na.rm = TRUE)
  } else {
    df$Immune_remodeling_score <- NA_real_
  }
  df
}

process_one_object <- function(row) {
  suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SummarizedExperiment)
  })
  decompressed_status <- "streamed_double_gzip_connection"
  obj <- read_double_gzip_rds(row$source_path)
  counts <- SummarizedExperiment::assay(obj, "counts")
  md <- as.data.frame(SingleCellExperiment::colData(obj))
  genes_upper <- toupper(rownames(obj))
  sc <- score_modules_by_cell(counts, genes_upper, sig_list)
  cell_scores <- bind_cols(
    tibble(
      dataset = "GSE266919",
      source_file = basename(row$source_path),
      source_celltype = row$source_celltype,
      cell_barcode = rownames(md)
    ),
    as_tibble(md),
    sc$scores
  )
  sample_scores <- cell_scores %>%
    group_by(dataset, source_celltype, Sample, Patient, Group, Tissue, Treatment, Efficacy, Response, MajorCluster) %>%
    summarise(
      n_cells = n(),
      across(all_of(names(sig_list)), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  coverage <- sc$coverage %>%
    mutate(
      dataset = "GSE266919",
      source_celltype = row$source_celltype,
      source_file = basename(row$source_path),
      decompression_status = decompressed_status,
      n_cells = ncol(counts),
      n_genes = nrow(counts),
      n_samples = dplyr::n_distinct(md$Sample),
      n_patients = dplyr::n_distinct(md$Patient),
      .before = 1
    )
  object_qc <- tibble(
    dataset = "GSE266919",
    source_celltype = row$source_celltype,
    source_file = basename(row$source_path),
    cache_path = "not_materialized_streamed_double_gzip",
    class = paste(class(obj), collapse = ";"),
    n_genes = nrow(counts),
    n_cells = ncol(counts),
    n_samples = dplyr::n_distinct(md$Sample),
    n_patients = dplyr::n_distinct(md$Patient),
    metadata_columns = paste(names(md), collapse = ";"),
    status = "parsed",
    message = "SingleCellExperiment counts and metadata parsed."
  )
  rm(obj, counts, md, cell_scores)
  gc()
  list(sample_scores = sample_scores, coverage = coverage, object_qc = object_qc)
}

wilcox_or_block <- function(df, score_col, contrast_label) {
  df <- df %>% filter(!is.na(.data[[score_col]]), !is.na(Response), Response %in% c("R", "NR"))
  n_r <- sum(df$Response == "R")
  n_nr <- sum(df$Response == "NR")
  if (n_r < 3 || n_nr < 3) {
    return(tibble(
      module = score_col, contrast = contrast_label, n_samples = nrow(df),
      n_patients = dplyr::n_distinct(df$Patient), n_R = n_r, n_NR = n_nr,
      effect = NA_real_, p_value = NA_real_, test = "not_fit",
      notes = "Response comparison requires at least 3 R and 3 NR samples."
    ))
  }
  wt <- suppressWarnings(wilcox.test(df[[score_col]] ~ df$Response, exact = FALSE))
  med_r <- median(df[[score_col]][df$Response == "R"], na.rm = TRUE)
  med_nr <- median(df[[score_col]][df$Response == "NR"], na.rm = TRUE)
  tibble(
    module = score_col, contrast = contrast_label, n_samples = nrow(df),
    n_patients = dplyr::n_distinct(df$Patient), n_R = n_r, n_NR = n_nr,
    effect = med_r - med_nr, p_value = wt$p.value, test = "wilcoxon_sample_level_R_minus_NR",
    notes = "Sample-level pseudobulk comparison; cells are not treated as independent response units."
  )
}

paired_change_or_block <- function(df, score_col, contrast_label) {
  wide <- df %>%
    filter(Group %in% c("Pre-treatment", "Post-treatment")) %>%
    group_by(Patient, Group) %>%
    summarise(score = mean(.data[[score_col]], na.rm = TRUE), n_cells = sum(n_cells), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = score)
  if (!all(c("Pre-treatment", "Post-treatment") %in% names(wide))) {
    return(tibble(
      module = score_col, contrast = contrast_label, n_pairs = 0L,
      effect = NA_real_, p_value = NA_real_, test = "not_fit",
      notes = "No paired pre/post samples available."
    ))
  }
  wide <- wide %>% filter(!is.na(`Pre-treatment`), !is.na(`Post-treatment`)) %>%
    mutate(diff_post_minus_pre = `Post-treatment` - `Pre-treatment`)
  if (nrow(wide) < 5) {
    return(tibble(
      module = score_col, contrast = contrast_label, n_pairs = nrow(wide),
      effect = median(wide$diff_post_minus_pre, na.rm = TRUE), p_value = NA_real_,
      test = "descriptive_only",
      notes = "Paired treatment-change test requires at least 5 patients."
    ))
  }
  wt <- suppressWarnings(wilcox.test(wide$diff_post_minus_pre, mu = 0, exact = FALSE))
  tibble(
    module = score_col, contrast = contrast_label, n_pairs = nrow(wide),
    effect = median(wide$diff_post_minus_pre, na.rm = TRUE), p_value = wt$p.value,
    test = "paired_wilcoxon_patient_level_post_minus_pre",
    notes = "Patient-level paired test; cells are not treated as independent response units."
  )
}

make_cell_source_figure <- function(celltype_summary, path) {
  plot_df <- celltype_summary %>%
    filter(module %in% c("APC_CXCL9_axis", "HLAII", "CXCL9_FOLR2_macrophage", "DC_HLA_antigen_presentation",
                         "Cytolytic", "Proliferation", "Immune_remodeling_score")) %>%
    group_by(module) %>%
    mutate(display_z = safe_z(mean_score)) %>%
    ungroup()
  p <- ggplot(plot_df, aes(x = module, y = source_celltype, fill = display_z)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient2(low = "#2B6CB0", mid = "white", high = "#B83280", midpoint = 0, na.value = "grey90") +
    labs(x = NULL, y = NULL, fill = "z") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
  ggsave(path, p, width = 7.4, height = 3.8, device = "pdf")
}

make_response_figure <- function(sample_scores, path) {
  plot_df <- sample_scores %>%
    filter(Group == "Pre-treatment", Response %in% c("R", "NR"), !is.na(APC_CXCL9_axis))
  if (!nrow(plot_df)) {
    pdf(path, width = 7, height = 4)
    plot.new()
    text(0.5, 0.5, "No pre-treatment response-labeled sample-level scRNA data.")
    dev.off()
    return(invisible(FALSE))
  }
  p <- ggplot(plot_df, aes(x = Response, y = APC_CXCL9_axis, fill = Response)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = 1.5, alpha = 0.75) +
    facet_wrap(~ source_celltype, scales = "free_y") +
    scale_fill_manual(values = c(NR = "#9CA3AF", R = "#2B6CB0")) +
    labs(x = NULL, y = "Sample-level APC/CXCL9/HLA-II axis") +
    theme_bw(base_size = 9) +
    theme(legend.position = "none", panel.grid.minor = element_blank())
  ggsave(path, p, width = 7.4, height = 4.2, device = "pdf")
}

write_blocker <- function(object_qc, package_qc, reason) {
  lines <- c(
    "# Public scRNA Response or Blocker Report",
    "",
    paste("Generated:", Sys.time()),
    "",
    paste("Status:", reason),
    "",
    "## Object QC",
    "",
    paste(capture.output(print(object_qc, n = Inf)), collapse = "\n"),
    "",
    "## Package QC",
    "",
    paste(capture.output(print(package_qc, n = Inf)), collapse = "\n"),
    "",
    "No cell-level response tests were run."
  )
  writeLines(lines, file.path(reports_dir, "public_scrna_response_or_blocker_report.md"))
}

main <- function() {
  package_before <- pkg_status("before_analysis")
  install_qc <- ensure_single_cell_experiment()
  package_after <- pkg_status("after_analysis_package_check")
  readr::write_tsv(bind_rows(package_before, package_after), file.path(results_dir, "public_scrna_analysis_package_qc_external_validation.tsv"))
  readr::write_tsv(install_qc, file.path(results_dir, "public_scrna_singlecellexperiment_install_qc_external_validation.tsv"))
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    write_blocker(tibble(), bind_rows(package_before, package_after), "SingleCellExperiment unavailable; cannot parse GSE266919 RDS objects.")
    stop("SingleCellExperiment unavailable.")
  }

  accession_dir <- file.path(public_geo_root, "GSE266919")
  cache_dir <- file.path(accession_dir, "decompressed_cache")
  files <- tibble(
    source_celltype = c("Bcell", "CD4Tcell", "CD8Tcell", "Myeloid", "NKcell"),
    filename = c("GSE266919_Bcell.rds.gz", "GSE266919_CD4Tcell.rds.gz", "GSE266919_CD8Tcell.rds.gz",
                 "GSE266919_Myeloid.rds.gz", "GSE266919_NKcell.rds.gz")
  ) %>%
    mutate(
      source_path = file.path(accession_dir, filename),
      cache_path = file.path(cache_dir, sub("\\.gz$", "", filename))
    )
  missing_files <- files %>% filter(!file.exists(source_path))
  if (nrow(missing_files)) {
    object_qc <- missing_files %>%
      transmute(dataset = "GSE266919", source_celltype, source_file = filename, cache_path,
                class = NA_character_, n_genes = NA_integer_, n_cells = NA_integer_,
                n_samples = NA_integer_, n_patients = NA_integer_, metadata_columns = NA_character_,
                status = "missing_source_file", message = "Run script 23 acquisition first.")
    write_blocker(object_qc, bind_rows(package_before, package_after), "Missing GSE266919 source files.")
    stop("Missing GSE266919 files.")
  }

  processed <- lapply(seq_len(nrow(files)), function(i) {
    row <- files[i, ]
    tryCatch(process_one_object(row), error = function(e) {
      list(
        sample_scores = tibble(),
        coverage = tibble(),
        object_qc = tibble(
          dataset = "GSE266919", source_celltype = row$source_celltype,
          source_file = basename(row$source_path), cache_path = row$cache_path,
          class = NA_character_, n_genes = NA_integer_, n_cells = NA_integer_,
          n_samples = NA_integer_, n_patients = NA_integer_, metadata_columns = NA_character_,
          status = "parse_failed", message = one_line(conditionMessage(e))
        )
      )
    })
  })

  object_qc <- bind_rows(lapply(processed, `[[`, "object_qc"))
  coverage <- bind_rows(lapply(processed, `[[`, "coverage"))
  sample_scores <- bind_rows(lapply(processed, `[[`, "sample_scores"))
  readr::write_tsv(object_qc, file.path(results_dir, "public_scrna_object_qc_external_validation.tsv"))
  readr::write_tsv(coverage, file.path(results_dir, "public_scrna_signature_coverage_external_validation.tsv"))

  if (!nrow(sample_scores)) {
    write_blocker(object_qc, bind_rows(package_before, package_after), "No GSE266919 object could be parsed into sample-level module scores.")
    readr::write_tsv(tibble(), file.path(results_dir, "public_scrna_response_module_by_celltype.tsv"))
    readr::write_tsv(tibble(), file.path(results_dir, "public_scrna_patient_level_module_tests.tsv"))
    stop("No parsed scRNA sample scores.")
  }

  sample_scores <- add_composites(sample_scores)
  score_cols <- c(names(sig_list), "APC_CXCL9_axis", "Immune_remodeling_score")
  readr::write_tsv(sample_scores, file.path(results_dir, "public_scrna_patient_level_pseudobulk_scores.tsv"))

  celltype_summary <- sample_scores %>%
    select(dataset, source_celltype, Sample, Patient, Group, Response, n_cells, all_of(score_cols)) %>%
    pivot_longer(cols = all_of(score_cols), names_to = "module", values_to = "score") %>%
    group_by(dataset, source_celltype, module) %>%
    summarise(
      n_samples = n_distinct(Sample),
      n_patients = n_distinct(Patient),
      n_cells = sum(n_cells, na.rm = TRUE),
      mean_score = mean(score, na.rm = TRUE),
      median_score = median(score, na.rm = TRUE),
      sd_score = sd(score, na.rm = TRUE),
      .groups = "drop"
    )
  readr::write_tsv(celltype_summary, file.path(results_dir, "public_scrna_response_module_by_celltype.tsv"))

  tests <- list()
  for (ct in unique(sample_scores$source_celltype)) {
    df_ct <- sample_scores %>% filter(source_celltype == ct)
    for (module in c("APC_CXCL9_axis", "Immune_remodeling_score", "HLAII", "CXCL9_FOLR2_macrophage",
                     "DC_HLA_antigen_presentation", "Cytolytic", "Proliferation")) {
      tests[[length(tests) + 1]] <- df_ct %>%
        filter(Group == "Pre-treatment") %>%
        wilcox_or_block(module, paste0(ct, ": pre-treatment R vs NR")) %>%
        mutate(dataset = "GSE266919", source_celltype = ct, comparison_level = "sample_level_pretreatment_response", .before = 1)
      tests[[length(tests) + 1]] <- df_ct %>%
        paired_change_or_block(module, paste0(ct, ": post-treatment minus pre-treatment")) %>%
        mutate(dataset = "GSE266919", source_celltype = ct, comparison_level = "patient_level_paired_treatment_change", .before = 1)
    }
  }
  tests_tbl <- bind_rows(tests)
  readr::write_tsv(tests_tbl, file.path(results_dir, "public_scrna_patient_level_module_tests.tsv"))

  make_cell_source_figure(celltype_summary, file.path(figures_dir, "public_scrna_cell_source_APC_CXCL9_axis.pdf"))
  make_response_figure(sample_scores, file.path(figures_dir, "public_scrna_response_APC_CXCL9_axis_pretreatment.pdf"))

  response_axis <- tests_tbl %>%
    filter(comparison_level == "sample_level_pretreatment_response", module == "APC_CXCL9_axis") %>%
    select(source_celltype, n_samples, n_patients, n_R, n_NR, effect, p_value, test, notes)
  source_axis <- celltype_summary %>%
    filter(module == "APC_CXCL9_axis") %>%
    arrange(desc(mean_score)) %>%
    select(source_celltype, n_samples, n_patients, n_cells, mean_score, median_score)
  lines <- c(
    "# Public scRNA Response or Blocker Report",
    "",
    paste("Generated:", Sys.time()),
    "",
    "GSE266919 was successfully acquired and parsed as double-gzipped SingleCellExperiment RDS objects. All response analyses were performed after aggregating module scores to sample/patient-level units; cells were not treated as independent response observations.",
    "",
    "## Cell-source APC/CXCL9/HLA-II Axis Summary",
    "",
    paste(capture.output(print(source_axis, n = Inf)), collapse = "\n"),
    "",
    "## Pretreatment Response Tests",
    "",
    paste(capture.output(print(response_axis, n = Inf)), collapse = "\n"),
    "",
    "## Limitations",
    "",
    "- These analyses support treatment-context cell-source biology, not an independent clinical prediction claim.",
    "- Paired pre/post tests are descriptive where patient-pair counts are low."
  )
  writeLines(lines, file.path(reports_dir, "public_scrna_response_or_blocker_report.md"))

  message("Round 4 public TNBC scRNA cell-source and response analysis complete.")
  invisible(TRUE)
}

main()
