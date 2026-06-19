#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(yaml)
})

root <- getwd()
cfg <- yaml::read_yaml(file.path(root, "config", "project_config.yaml"))
sig_path <- file.path(root, cfg$paths$signatures_csv)
if (!file.exists(sig_path)) stop("Signature file not found: ", sig_path)

sig <- readr::read_csv(sig_path, show_col_types = FALSE) %>%
  mutate(
    gene = toupper(gene),
    signature = as.character(signature),
    module_group = as.character(module_group)
  ) %>%
  distinct(signature, gene, .keep_all = TRUE)

# Basic sanity checks
summary_tbl <- sig %>%
  group_by(signature, module_group, role) %>%
  summarise(n_genes = n_distinct(gene), genes = paste(sort(unique(gene)), collapse = ";"), .groups = "drop") %>%
  arrange(module_group, signature)

# Identify genes duplicated across signatures. This is allowed for HLA and controls but should be reported.
dup_tbl <- sig %>%
  count(gene, sort = TRUE) %>%
  filter(n > 1) %>%
  left_join(sig %>% group_by(gene) %>% summarise(signatures = paste(unique(signature), collapse = ";"), .groups = "drop"), by = "gene")

# Empty coverage template for downstream completion by dataset.
coverage_template <- summary_tbl %>%
  transmute(dataset = NA_character_, signature, n_signature_genes = n_genes,
            n_genes_found = NA_integer_, coverage = NA_real_, missing_genes = NA_character_)

out_clean <- file.path(root, "data", "processed", "TNBC_BIR_signatures_clean.tsv")
out_summary <- file.path(root, "results", "signature_summary.tsv")
out_dup <- file.path(root, "results", "signature_duplicate_genes.tsv")
out_cov <- file.path(root, "results", "signature_gene_coverage_template.tsv")

readr::write_tsv(sig, out_clean)
readr::write_tsv(summary_tbl, out_summary)
readr::write_tsv(dup_tbl, out_dup)
readr::write_tsv(coverage_template, out_cov)

message("Prepared signatures:")
print(summary_tbl)
message("Wrote: ", out_clean)
