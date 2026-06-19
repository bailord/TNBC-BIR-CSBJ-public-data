# Utility functions for TNBC-BIR signature scoring.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

make_unique_gene_matrix <- function(expr, duplicate_rule = c("highest_mad", "mean")) {
  duplicate_rule <- match.arg(duplicate_rule)
  stopifnot(!is.null(rownames(expr)))
  rownames(expr) <- toupper(rownames(expr))
  expr <- as.matrix(expr)
  expr <- expr[!is.na(rownames(expr)) & rownames(expr) != "", , drop = FALSE]

  if (!any(duplicated(rownames(expr)))) return(expr)

  if (duplicate_rule == "mean") {
    df <- as.data.frame(expr) %>% rownames_to_column("gene")
    out <- df %>% group_by(gene) %>% summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
    mat <- as.matrix(out[, -1])
    rownames(mat) <- out$gene
    return(mat)
  }

  mad_tbl <- tibble(gene = rownames(expr), idx = seq_len(nrow(expr)), mad = apply(expr, 1, mad, na.rm = TRUE)) %>%
    group_by(gene) %>%
    slice_max(mad, n = 1, with_ties = FALSE) %>%
    ungroup()
  expr[mad_tbl$idx, , drop = FALSE]
}

maybe_log2 <- function(expr) {
  expr <- as.matrix(expr)
  q <- stats::quantile(expr, probs = c(0.5, 0.9, 0.99), na.rm = TRUE)
  # Heuristic: raw counts/intensities are often large; log-scale microarray/RNAseq values are usually < 30.
  if (q[[3]] > 100 || q[[2]] > 50) {
    message("Applying log2(x + 1) transform based on expression quantiles.")
    expr <- log2(expr + 1)
  }
  expr
}

zscore_rows <- function(expr) {
  z <- t(scale(t(expr)))
  z[is.na(z)] <- 0
  z
}

score_signatures <- function(expr, sig_tbl, min_coverage = 0.5, duplicate_rule = "highest_mad") {
  expr <- make_unique_gene_matrix(expr, duplicate_rule = duplicate_rule)
  expr <- maybe_log2(expr)
  expr_z <- zscore_rows(expr)

  sig_list <- split(toupper(sig_tbl$gene), sig_tbl$signature)
  scores <- list()
  coverage <- list()

  for (nm in names(sig_list)) {
    genes <- unique(sig_list[[nm]])
    present <- intersect(genes, rownames(expr_z))
    cov <- length(present) / length(genes)
    missing <- setdiff(genes, rownames(expr_z))
    if (cov >= min_coverage && length(present) > 0) {
      scores[[nm]] <- colMeans(expr_z[present, , drop = FALSE], na.rm = TRUE)
    } else {
      scores[[nm]] <- rep(NA_real_, ncol(expr_z))
      names(scores[[nm]]) <- colnames(expr_z)
    }
    coverage[[nm]] <- tibble(
      signature = nm,
      n_signature_genes = length(genes),
      n_genes_found = length(present),
      coverage = cov,
      missing_genes = paste(missing, collapse = ";")
    )
  }

  score_df <- bind_rows(lapply(names(scores), function(nm) tibble(sample_id = names(scores[[nm]]), signature = nm, score = as.numeric(scores[[nm]])))) %>%
    pivot_wider(names_from = signature, values_from = score)

  # Composite score definitions. Use scale() on module scores across samples.
  zcol <- function(x) as.numeric(scale(x))
  if (all(c("Mast_cell", "Tfh_TLS_like", "CD8_Tsem_memory", "CXCL9_FOLR2_macrophage", "DC_HLA_antigen_presentation") %in% colnames(score_df))) {
    score_df <- score_df %>%
      mutate(
        Immune_remodeling_score = rowMeans(across(c(Mast_cell, Tfh_TLS_like, CD8_Tsem_memory, CXCL9_FOLR2_macrophage, DC_HLA_antigen_presentation)), na.rm = TRUE)
      )
  }
  if (all(c("Immune_remodeling_score", "Basoluminal_plasticity", "Suppressive_myeloid") %in% colnames(score_df))) {
    score_df <- score_df %>%
      mutate(
        Basoluminal_score = Basoluminal_plasticity,
        Suppressive_myeloid_score = Suppressive_myeloid,
        TNBC_BIR_favorable_score = zcol(Immune_remodeling_score) - zcol(Basoluminal_score) - zcol(Suppressive_myeloid_score)
      )
  }

  coverage_df <- bind_rows(coverage)
  list(scores = score_df, coverage = coverage_df)
}

assign_bir_quadrants <- function(score_df) {
  stopifnot(all(c("Basoluminal_score", "Immune_remodeling_score") %in% colnames(score_df)))
  score_df %>%
    mutate(
      Basoluminal_high = Basoluminal_score >= median(Basoluminal_score, na.rm = TRUE),
      Immune_high = Immune_remodeling_score >= median(Immune_remodeling_score, na.rm = TRUE),
      Suppressive_high = ifelse("Suppressive_myeloid_score" %in% names(.), Suppressive_myeloid_score >= median(Suppressive_myeloid_score, na.rm = TRUE), NA),
      BIR_quadrant = case_when(
        Immune_high & !Basoluminal_high ~ "Immune-remodeled",
        !Immune_high & Basoluminal_high ~ "Basoluminal-cold",
        Immune_high & Basoluminal_high ~ "Basoluminal-inflamed",
        TRUE ~ "Immune-cold/non-basoluminal"
      )
    )
}
