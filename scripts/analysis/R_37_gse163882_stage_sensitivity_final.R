#!/usr/bin/env Rscript

set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("reports", showWarnings = FALSE, recursive = TRUE)

scores_path <- "data/processed/external_validation_external_NAC_scores_clinical.tsv"
series_path <- "<PUBLIC_BIOINF_DATABASES>/GEO/GSE163882/GSE163882_series_matrix.txt.gz"

safe_write_blocker <- function(reason) {
  out <- tibble(
    dataset = "GSE163882",
    analysis_subset = "primary_TNBC",
    endpoint = "pCR",
    n_total = NA_real_,
    n_event_or_pCR = NA_real_,
    model_name = "stage_sensitivity_blocker",
    variable = NA_character_,
    effect_type = NA_character_,
    estimate = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    p_value = NA_real_,
    covariates = NA_character_,
    notes = reason
  )
  write_tsv(out, "results/GSE163882_stage_sensitivity_models.tsv")
  writeLines(c(
    "# GSE163882 stage sensitivity report",
    "",
    "Stage sensitivity could not be completed.",
    "",
    paste0("Blocker: ", reason),
    "",
    "Round 5 interpretation is retained: supportive but proliferation-sensitive evidence."
  ), "reports/GSE163882_stage_sensitivity_report.md")
  stop(reason, call. = FALSE)
}

if (!file.exists(scores_path)) {
  safe_write_blocker(paste("Missing score/clinical table:", scores_path))
}
if (!file.exists(series_path)) {
  safe_write_blocker(paste("Missing GEO series matrix:", series_path))
}

stripq <- function(x) gsub('^"|"$', "", x)

lines <- readLines(gzfile(series_path), warn = FALSE)
geo_line <- lines[str_detect(lines, "^!Sample_geo_accession\\t")]
title_line <- lines[str_detect(lines, "^!Sample_title\\t")]
stage_line <- lines[str_detect(lines, "^!Sample_characteristics_ch1\\t") & str_detect(lines, "breast cancer stage:")]

if (length(geo_line) != 1 || length(title_line) != 1 || length(stage_line) != 1) {
  safe_write_blocker("Could not uniquely parse !Sample_geo_accession, !Sample_title, and breast cancer stage rows from series matrix.")
}

geo_ids <- stripq(strsplit(geo_line, "\t", fixed = TRUE)[[1]][-1])
sample_titles <- stripq(strsplit(title_line, "\t", fixed = TRUE)[[1]][-1])
stage_values <- stripq(strsplit(stage_line, "\t", fixed = TRUE)[[1]][-1])

if (length(geo_ids) != length(stage_values) || length(sample_titles) != length(stage_values)) {
  safe_write_blocker(paste("Stage row length mismatch:", length(geo_ids), "GSM IDs,", length(sample_titles), "titles vs", length(stage_values), "stage values."))
}

stage_tbl <- tibble(
  geo_sample_id = geo_ids,
  title = sample_titles,
  sample_id = str_replace(sample_titles, ":.*$", ""),
  stage_raw = str_replace(stage_values, "^breast cancer stage:\\s*", ""),
  stage_num = suppressWarnings(as.numeric(stage_raw)),
  stage_missing = is.na(stage_num),
  stage_iv = !is.na(stage_num) & stage_num == 4,
  stage_group = case_when(
    is.na(stage_num) ~ NA_character_,
    stage_num == 4 ~ "IV",
    TRUE ~ paste0("I-III_stage_", stage_num)
  )
)

scores <- read_tsv(scores_path, show_col_types = FALSE, progress = FALSE) %>%
  filter(dataset == "GSE163882", primary_tnbc %in% TRUE, pCR_mapped %in% TRUE) %>%
  mutate(
    pCR_binary = as.integer(pCR_binary),
    APC_CXCL9_axis = as.numeric(APC_CXCL9_axis),
    Proliferation = as.numeric(Proliferation),
    TNBC_BIR_favorable_score = as.numeric(TNBC_BIR_favorable_score)
  )

joined <- scores %>%
  left_join(
    stage_tbl %>% select(stage_geo_sample_id = geo_sample_id, stage_title = title, sample_id, stage_raw, stage_num, stage_missing, stage_iv, stage_group),
    by = "sample_id"
  )

join_n <- sum(!is.na(joined$stage_raw))
if (nrow(joined) == 0 || join_n < 80) {
  safe_write_blocker(paste0("Stage join insufficient: joined ", join_n, " of ", nrow(joined), " primary TNBC samples."))
}

fit_model <- function(df, model_name, formula, variables, notes) {
  cc <- df %>%
    filter(!is.na(pCR_binary), !is.na(APC_CXCL9_axis)) %>%
    mutate(stage_group = factor(stage_group))
  if (nrow(cc) < 25 || sum(cc$pCR_binary == 1) < 8 || sum(cc$pCR_binary == 0) < 8) {
    return(tibble(
      dataset = "GSE163882",
      analysis_subset = "primary_TNBC",
      endpoint = "pCR",
      n_total = nrow(cc),
      n_event_or_pCR = sum(cc$pCR_binary == 1),
      model_name = model_name,
      variable = variables,
      effect_type = "OR",
      estimate = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      covariates = deparse(formula),
      notes = paste(notes, "Not fit because subset is underpowered.")
    ))
  }
  fit <- tryCatch(glm(formula, data = cc, family = binomial()), error = function(e) e)
  if (inherits(fit, "error")) {
    return(tibble(
      dataset = "GSE163882",
      analysis_subset = "primary_TNBC",
      endpoint = "pCR",
      n_total = nrow(cc),
      n_event_or_pCR = sum(cc$pCR_binary == 1),
      model_name = model_name,
      variable = variables,
      effect_type = "OR",
      estimate = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      covariates = deparse(formula),
      notes = paste(notes, "Model failed:", fit$message)
    ))
  }
  co <- summary(fit)$coefficients
  out <- lapply(variables, function(v) {
    coef_name <- rownames(co)[rownames(co) == v]
    if (length(coef_name) == 0) coef_name <- rownames(co)[str_detect(rownames(co), fixed(v))]
    if (length(coef_name) == 0) {
      return(tibble(variable = v, estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_))
    }
    beta <- co[coef_name[1], "Estimate"]
    se <- co[coef_name[1], "Std. Error"]
    tibble(
      variable = coef_name[1],
      estimate = exp(beta),
      conf_low = exp(beta - 1.96 * se),
      conf_high = exp(beta + 1.96 * se),
      p_value = co[coef_name[1], "Pr(>|z|)"]
    )
  }) %>% bind_rows()
  out %>%
    mutate(
      dataset = "GSE163882",
      analysis_subset = "primary_TNBC",
      endpoint = "pCR",
      n_total = nrow(model.frame(fit)),
      n_event_or_pCR = sum(model.frame(fit)$pCR_binary == 1),
      model_name = model_name,
      effect_type = "OR",
      covariates = deparse(formula),
      notes = notes
    ) %>%
    select(dataset, analysis_subset, endpoint, n_total, n_event_or_pCR, model_name,
           variable, effect_type, estimate, conf_low, conf_high, p_value, covariates, notes)
}

model_rows <- bind_rows(
  fit_model(
    joined,
    "all_stages_axis_only",
    pCR_binary ~ scale(APC_CXCL9_axis),
    "scale(APC_CXCL9_axis)",
    "All primary TNBC samples; stage parsed but not filtered."
  ),
  fit_model(
    joined %>% filter(!stage_iv),
    "exclude_stage_IV_axis_only",
    pCR_binary ~ scale(APC_CXCL9_axis),
    "scale(APC_CXCL9_axis)",
    "Stage IV samples excluded; missing stage retained."
  ),
  fit_model(
    joined %>% filter(!stage_missing),
    "exclude_missing_stage_axis_only",
    pCR_binary ~ scale(APC_CXCL9_axis),
    "scale(APC_CXCL9_axis)",
    "Missing stage excluded; stage IV retained."
  ),
  fit_model(
    joined %>% filter(!stage_missing, !stage_iv),
    "exclude_stage_IV_and_missing_axis_only",
    pCR_binary ~ scale(APC_CXCL9_axis),
    "scale(APC_CXCL9_axis)",
    "Stage IV and missing-stage samples excluded."
  ),
  fit_model(
    joined %>% filter(!stage_missing),
    "stage_adjusted_axis_plus_proliferation",
    pCR_binary ~ scale(APC_CXCL9_axis) + scale(Proliferation) + stage_group,
    c("scale(APC_CXCL9_axis)", "scale(Proliferation)"),
    "Exploratory complete-case stage-adjusted model; interpret as sensitivity only."
  )
)

write_tsv(model_rows, "results/GSE163882_stage_sensitivity_models.tsv")

qc_tbl <- joined %>%
  summarise(
    dataset = "GSE163882",
    all_primary_tnbc_pcr_mapped = n(),
    stage_joined = sum(!is.na(stage_raw)),
    stage_missing = sum(stage_missing %in% TRUE),
    stage_IV = sum(stage_iv %in% TRUE),
    pCR = sum(pCR_binary == 1, na.rm = TRUE),
    RD = sum(pCR_binary == 0, na.rm = TRUE)
  )

stage_dist <- joined %>%
  mutate(stage_display = if_else(is.na(stage_raw), "unjoined", stage_raw)) %>%
  count(stage_display, pCR_binary, name = "n") %>%
  arrange(stage_display, pCR_binary)

axis_all <- model_rows %>%
  filter(model_name == "all_stages_axis_only", variable == "scale(APC_CXCL9_axis)") %>%
  slice_head(n = 1)
axis_adj <- model_rows %>%
  filter(model_name == "stage_adjusted_axis_plus_proliferation", variable == "scale(APC_CXCL9_axis)") %>%
  slice_head(n = 1)

interpretation <- "supportive but proliferation-sensitive"
if (nrow(axis_all) > 0 && nrow(axis_adj) > 0 &&
    !is.na(axis_all$p_value[1]) && !is.na(axis_adj$p_value[1]) &&
    axis_all$p_value[1] < 0.05 && axis_adj$p_value[1] >= 0.05) {
  interpretation <- "supportive but proliferation-sensitive; stage sensitivity does not promote this to proliferation-independent evidence"
}
if (nrow(axis_all) > 0 && !is.na(axis_all$estimate[1]) && axis_all$estimate[1] < 1) {
  interpretation <- "exploratory supplement because stage sensitivity reverses the axis direction"
}

report <- c(
  "# GSE163882 stage sensitivity report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Join QC",
  "",
  paste(capture.output(print(as.data.frame(qc_tbl), row.names = FALSE)), collapse = "\n"),
  "",
  "## Stage distribution by pCR status",
  "",
  paste(capture.output(print(as.data.frame(stage_dist), row.names = FALSE)), collapse = "\n"),
  "",
  "## Interpretation lock",
  "",
  paste0("- Final interpretation: ", interpretation, "."),
  "- Stage is parsed from the GEO series matrix `breast cancer stage` row and joined by `geo_sample_id`.",
  "- Stage-adjusted models are complete-case exploratory sensitivity analyses, not upgraded biomarker claims.",
  "",
  "## Output",
  "",
  "- `results/GSE163882_stage_sensitivity_models.tsv`"
)

writeLines(report, "reports/GSE163882_stage_sensitivity_report.md")

cat("Wrote GSE163882 stage sensitivity with", nrow(joined), "primary TNBC samples; stage joined", join_n, "\n")
