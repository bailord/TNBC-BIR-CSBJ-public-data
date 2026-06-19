#!/usr/bin/env Rscript

set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

dir.create("figures/final", showWarnings = FALSE, recursive = TRUE)

read_tsv_if <- function(path) {
  if (!file.exists(path)) return(tibble())
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

save_plot <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = "pdf")
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 300)
}

wong <- c(
  orange = "#E69F00",
  sky = "#56B4E9",
  green = "#009E73",
  blue = "#0072B2",
  vermillion = "#D55E00",
  purple = "#CC79A7",
  black = "#000000",
  grey = "#7A7A7A"
)

theme_pub <- function(base_size = 8) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "sans"),
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

table2 <- read_tsv_if("results/final_table2_primary_cross_cohort_effects.tsv")
table3 <- read_tsv_if("results/final_table3_scRNA_source_summary.tsv")

# Figure 1: locked framework
nodes <- tibble(
  id = c("Locked signatures", "Tumor axes", "Immune-reactive axis", "Controls", "Context tests", "Claim boundary"),
  x = c(1, 2.4, 2.4, 2.4, 3.8, 5.2),
  y = c(2, 3, 2, 1, 2, 2),
  label = c(
    "Locked TNBC-BIR\nsignatures",
    "Basoluminal\nplasticity",
    "APC/CXCL9/HLA-II\nimmune-reactive axis",
    "Proliferation,\ncytolytic controls",
    "I-SPY2 pCR,\nexternal NAC,\nMETABRIC,\nscRNA source",
    "State decomposition;\nnot a universal\npCR predictor"
  ),
  fill = c("#F4F4F4", "#FDEBD0", "#D6EAF8", "#E8F6EF", "#EBDEF0", "#F4F4F4")
)
edges <- tibble(
  x = c(1.45, 1.45, 1.45, 2.85, 4.25),
  y = c(2, 2, 2, 2, 2),
  xend = c(2.0, 2.0, 2.0, 3.35, 4.75),
  yend = c(3, 2, 1, 2, 2)
)
fig1 <- ggplot() +
  geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.12, "inches")), linewidth = 0.4, color = "grey35") +
  geom_rect(data = nodes, aes(xmin = x - 0.45, xmax = x + 0.45, ymin = y - 0.35, ymax = y + 0.35, fill = fill),
            color = "grey30", linewidth = 0.35) +
  geom_text(data = nodes, aes(x = x, y = y, label = label), size = 2.8, lineheight = 0.95) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0.45, 5.75), ylim = c(0.45, 3.55), expand = FALSE) +
  labs(title = "Locked TNBC state-decomposition framework") +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11, margin = margin(b = 8)))
save_plot(fig1, "figures/final/Figure1_locked_TNBC_BIR_framework", 7.0, 3.7)

# Figure 2: manuscript-facing forest, one row per dataset/context
key_rows <- table2 %>%
  filter(score_or_axis == "APC_CXCL9_axis") %>%
  mutate(
    row_key = case_when(
      dataset == "GSE194040" & model_name == "axis_treatment_adjusted" ~ "GSE194040 I-SPY2 pCR",
      dataset == "GSE25066_RMA" & endpoint == "pCR" & model_name == "pCR_APC_CXCL9_axis" ~ "GSE25066 RMA pCR",
      dataset == "GSE25066_RMA" & endpoint == "DRFS" ~ "GSE25066 RMA DRFS",
      dataset == "GSE163882" & source_table == "GSE163882_axis_models_final" & model_name == "axis_only" ~ "GSE163882 pCR",
      dataset == "GSE41998" & model_name == "axis_only" ~ "GSE41998 pCR",
      dataset == "GSE20194" & model_name == "axis_only" ~ "GSE20194 pCR",
      dataset == "manual_external_NAC_meta" ~ "Manual external NAC meta",
      dataset == "METABRIC" & model_name == "axis_age_adjusted" ~ "METABRIC RFS",
      TRUE ~ NA_character_
    ),
    context_class = case_when(
      dataset == "GSE194040" ~ "I-SPY2 context",
      dataset == "GSE25066_RMA" ~ "valid negative",
      dataset == "GSE163882" ~ "supportive/sensitive",
      dataset %in% c("GSE41998", "GSE20194", "manual_external_NAC_meta") ~ "external directional",
      dataset == "METABRIC" ~ "recurrence biology",
      TRUE ~ "context"
    )
  ) %>%
  filter(!is.na(row_key), !is.na(estimate), !is.na(conf_low), !is.na(conf_high)) %>%
  distinct(row_key, .keep_all = TRUE) %>%
  mutate(
    row_key = factor(row_key, levels = rev(c(
      "GSE194040 I-SPY2 pCR", "GSE25066 RMA pCR", "GSE25066 RMA DRFS",
      "GSE163882 pCR", "GSE41998 pCR", "GSE20194 pCR",
      "Manual external NAC meta", "METABRIC RFS"
    ))),
    effect_label = if_else(effect_type == "HR", "HR per 1 SD", "OR per 1 SD")
  )

fig2 <- key_rows %>%
  ggplot(aes(x = estimate, y = row_key, color = context_class)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.18, linewidth = 0.55) +
  geom_point(size = 2.2) +
  scale_x_log10() +
  scale_color_manual(values = c(
    "I-SPY2 context" = unname(wong["blue"]),
    "valid negative" = unname(wong["black"]),
    "supportive/sensitive" = unname(wong["vermillion"]),
    "external directional" = unname(wong["green"]),
    "recurrence biology" = unname(wong["purple"]),
    "context" = unname(wong["grey"])
  ), name = NULL) +
  labs(
    x = "Effect estimate (log scale)",
    y = NULL,
    title = "Cross-cohort APC/CXCL9/HLA-II immune-reactive axis effects",
    subtitle = "One manuscript-facing row per dataset/context"
  ) +
  theme_pub(8)
save_plot(fig2, "figures/final/Figure2_cross_cohort_APC_CXCL9_effects", 7.0, 4.8)

# Figure 3: context map without overplotted p-values
context_map <- table2 %>%
  filter(score_or_axis %in% c("APC_CXCL9_axis", "Proliferation", "TNBC_BIR_favorable_score_descriptive")) %>%
  mutate(
    dataset_context = case_when(
      dataset == "GSE25066_RMA" ~ paste(dataset, endpoint, sep = " "),
      dataset == "METABRIC" ~ "METABRIC RFS",
      dataset == "manual_external_NAC_meta" ~ "external NAC meta",
      TRUE ~ dataset
    ),
    direction = case_when(
      is.na(estimate) ~ "not estimated",
      effect_type == "HR" & estimate < 1 ~ "favorable direction",
      effect_type == "OR" & estimate > 1 ~ "higher pCR odds",
      effect_type == "OR" & estimate <= 1 ~ "null/opposite",
      effect_type == "HR" & estimate >= 1 ~ "null/opposite",
      TRUE ~ "context"
    ),
    signal_strength = case_when(
      !is.na(p_value) & p_value < 0.05 ~ "nominal p<0.05",
      TRUE ~ "not nominal"
    )
  ) %>%
  group_by(dataset_context, score_or_axis) %>%
  slice_min(order_by = p_value, n = 1, with_ties = FALSE) %>%
  ungroup()

fig3 <- context_map %>%
  ggplot(aes(x = score_or_axis, y = dataset_context, fill = direction, shape = signal_strength)) +
  geom_tile(color = "white", linewidth = 0.35, alpha = 0.85) +
  geom_point(size = 2.2, color = "grey20") +
  scale_fill_manual(values = c(
    "higher pCR odds" = "#D6EAF8",
    "favorable direction" = "#D5F5E3",
    "null/opposite" = "#F4F4F4",
    "not estimated" = "#FFFFFF",
    "context" = "#F4F4F4"
  ), name = "Direction") +
  scale_shape_manual(values = c("nominal p<0.05" = 16, "not nominal" = 1), name = "Evidence") +
  labs(
    x = NULL,
    y = NULL,
    title = "Context-dependent axis and proliferation map",
    subtitle = "P-values are encoded by marker style to avoid label crowding"
  ) +
  theme_pub(8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(fig3, "figures/final/Figure3_axis_vs_proliferation_context_map", 7.0, 4.6)

# Figure 4: combined scRNA source and dynamics
major_source <- read_tsv_if("results/public_scrna_response_module_by_celltype.tsv") %>%
  filter(dataset == "GSE266919",
         source_celltype %in% c("Bcell", "Myeloid", "CD4Tcell", "CD8Tcell", "NKcell"),
         module %in% c("HLAII", "CXCL9_FOLR2_macrophage", "DC_HLA_antigen_presentation", "APC_CXCL9_axis")) %>%
  mutate(
    source_celltype = factor(source_celltype, levels = c("Bcell", "Myeloid", "CD4Tcell", "CD8Tcell", "NKcell")),
    module = factor(module, levels = c("HLAII", "CXCL9_FOLR2_macrophage", "DC_HLA_antigen_presentation", "APC_CXCL9_axis"))
  )
sub_source <- read_tsv_if("results/GSE266919_final_module_specific_source.tsv") %>%
  filter(source_level == "selected_subcluster",
         module %in% c("HLAII", "CXCL9_FOLR2_macrophage", "DC_HLA_antigen_presentation", "APC_CXCL9_axis"),
         source_celltype %in% c("Myeloid", "Bcell")) %>%
  group_by(source) %>%
  slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(source_celltype, desc(score)) %>%
  slice_head(n = 12)
dyn <- read_tsv_if("results/GSE266919_final_selected_subcluster_fraction_dynamics.tsv") %>%
  mutate(Subset = factor(Subset, levels = rev(unique(Subset))))

p4a <- major_source %>%
  ggplot(aes(x = module, y = source_celltype, fill = mean_score)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#56B4E9", mid = "white", high = "#D55E00", midpoint = 0, name = "Mean") +
  labs(x = NULL, y = NULL, title = "A. Major cell types") +
  theme_pub(7) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")

p4b <- sub_source %>%
  mutate(source = factor(source, levels = rev(source))) %>%
  ggplot(aes(x = score, y = source, fill = source_celltype)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(Myeloid = unname(wong["blue"]), Bcell = unname(wong["orange"])), name = NULL) +
  labs(x = "Top module score", y = NULL, title = "B. Selected APC/myeloid/B-cell subclusters") +
  theme_pub(7)

p4c <- dyn %>%
  ggplot(aes(x = median_delta, y = Subset)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(aes(color = p_adj_BH < 0.05), size = 1.8) +
  geom_text(aes(label = paste0("n=", n_pairs)), hjust = -0.25, size = 2.0) +
  scale_color_manual(values = c(`TRUE` = unname(wong["vermillion"]), `FALSE` = unname(wong["blue"])), labels = c("not FDR-significant", "FDR < 0.05"), name = NULL) +
  coord_cartesian(clip = "off") +
  labs(x = "Median post-minus-pre fraction", y = NULL, title = "C. Paired fraction dynamics") +
  theme_pub(7) +
  theme(plot.margin = margin(5.5, 28, 5.5, 5.5))

fig4 <- (p4a / p4b / p4c) + plot_layout(heights = c(0.8, 1.0, 1.0))
save_plot(fig4, "figures/final/Figure4_GSE266919_scRNA_source_and_dynamics", 7.0, 8.8)

# Figure 5: METABRIC recurrence biology support
metabric_rows <- table2 %>%
  filter(dataset == "METABRIC", score_or_axis == "APC_CXCL9_axis", endpoint == "RFS", !is.na(estimate)) %>%
  mutate(
    model_label = case_when(
      model_name == "axis_age_adjusted" ~ "Age-adjusted",
      model_name == "axis_age_stage_adjusted" ~ "Age + stage complete-case",
      TRUE ~ model_name
    ),
    model_label = factor(model_label, levels = rev(unique(model_label)))
  )

if (nrow(metabric_rows) == 0) {
  metabric_rows <- tibble(
    model_label = factor("Age-adjusted"),
    estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
    p_value = NA_real_, notes = "METABRIC rows unavailable"
  )
}

fig5 <- metabric_rows %>%
  ggplot(aes(x = estimate, y = model_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.18, linewidth = 0.6, na.rm = TRUE) +
  geom_point(size = 2.4, color = wong["purple"], na.rm = TRUE) +
  scale_x_log10() +
  labs(
    x = "Hazard ratio per 1 SD APC/CXCL9 axis",
    y = NULL,
    title = "METABRIC recurrence biology support",
    subtitle = "Shown as recurrence biology, not stage-independent biomarker proof"
  ) +
  theme_pub(8)
save_plot(fig5, "figures/final/Figure5_METABRIC_recurrence_biology", 5.0, 2.6)

# Supplementary aggregate BIR null figure
bir_rows <- table2 %>%
  filter(score_or_axis == "TNBC_BIR_favorable_score_descriptive", !is.na(estimate))
if (nrow(bir_rows) > 0) {
  supp_bir <- bir_rows %>%
    mutate(label = paste(dataset, endpoint, model_name, sep = " | "),
           label = factor(label, levels = rev(unique(label)))) %>%
    ggplot(aes(x = estimate, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
    geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.18, linewidth = 0.5) +
    geom_point(size = 2, color = wong["grey"]) +
    scale_x_log10() +
    labs(x = "Effect estimate", y = NULL, title = "Supplementary aggregate BIR descriptive/null summary") +
    theme_pub(8)
  save_plot(supp_bir, "figures/final/Supplementary_aggregate_BIR_descriptive_null", 6.0, 2.8)
}

cat("Wrote final manuscript figures to figures/final\n")
