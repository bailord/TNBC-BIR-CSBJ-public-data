#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: render_v2_5_submission_figures.R <output_dir>")
out <- normalizePath(args[[1]], mustWork = TRUE)
src <- file.path(out, "figures", "source_data")
fig <- file.path(out, "figures", "main")
sup <- file.path(out, "figures", "supplementary")
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
dir.create(sup, recursive = TRUE, showWarnings = FALSE)

navy <- "#173B57"
blue <- "#2878A5"
teal <- "#2A8C82"
orange <- "#D97925"
light_blue <- "#DDEBF3"
light_teal <- "#DCEFEA"
light_orange <- "#F7E5D2"
light_gray <- "#EFF1F2"
dark_gray <- "#4D555B"
mid_gray <- "#98A0A6"

theme_pub <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.5, color = navy, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.3, color = dark_gray),
      axis.title = element_text(face = "bold", color = dark_gray),
      axis.text = element_text(color = dark_gray),
      panel.grid = element_blank(),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 10, 8, 10)
    )
}

save_all <- function(p, stem, width, height, dpi = 320, supplementary = FALSE) {
  d <- if (supplementary) sup else fig
  ggsave(file.path(d, paste0(stem, ".pdf")), p, width = width, height = height,
         units = "in", device = "pdf", bg = "white", useDingbats = FALSE)
  ggsave(file.path(d, paste0(stem, ".png")), p, width = width, height = height,
         units = "in", dpi = dpi, bg = "white")
  if (!requireNamespace("svglite", quietly = TRUE)) {
    stop("Package 'svglite' is required for submission-grade SVG export")
  }
  ggsave(file.path(d, paste0(stem, ".svg")), p, width = width, height = height,
         units = "in", device = svglite::svglite, bg = "white")
}

# Figure 1: experimental and technical design. The visual is a single evidence
# trajectory, not four decorative cards.
f1 <- read_tsv(file.path(src, "Figure1_design_source_data.tsv"), show_col_types = FALSE)
stage_x <- c(0.75, 3.05, 5.75, 9.15)

p1 <- ggplot() +
  annotate("segment", x = 0.7, xend = 11.55, y = 5.22, yend = 5.22,
           arrow = arrow(length = unit(0.11, "in"), type = "closed"),
           color = "#A7B0B6", linewidth = 0.65) +
  annotate("segment", x = c(2.65, 5.35, 8.85), xend = c(2.65, 5.35, 8.85),
           y = 0.85, yend = 5.85, color = "#D7DCDF", linewidth = 0.45) +
  annotate("point", x = stage_x, y = 5.22, shape = 21, size = 5.1,
           fill = c(navy, teal, blue, orange), color = "white", stroke = 0.85) +
  annotate("text", x = stage_x, y = 5.22, label = 1:4, color = "white",
           fontface = "bold", size = 3.4) +
  annotate("text", x = c(1.42, 4.0, 7.08, 10.25), y = 5.72,
           label = c("Literature-anchored modules", "Cohort-wise fixed scoring",
                     "Evidence hierarchy", "Interpretation boundary"),
           color = navy, fontface = "bold", size = 4.05) +

  annotate("text", x = 0.32, y = 4.62, label = "Initial state decomposition",
           hjust = 0, color = dark_gray, fontface = "bold", size = 3.15) +
  annotate("segment", x = 0.34, xend = 0.34, y = 2.82, yend = 4.28,
           color = blue, linewidth = 2.4, lineend = "round") +
  annotate("text", x = 0.5, y = 4.17,
           label = "Basoluminal / epithelial\nImmune remodeling\nSuppressive myeloid\nProliferation / cytolytic",
           hjust = 0, vjust = 1, color = dark_gray, size = 3.0, lineheight = 1.25) +
  annotate("text", x = 0.32, y = 2.32, label = "Antigen-presentation summary",
           hjust = 0, color = navy, fontface = "bold", size = 3.15) +
  annotate("segment", x = 0.34, xend = 0.34, y = 1.02, yend = 2.03,
           color = teal, linewidth = 2.4, lineend = "round") +
  annotate("text", x = 0.5, y = 1.92,
           label = "HLA-II\nCXCL9 / FOLR2 macrophage\nDendritic-cell / HLA",
           hjust = 0, vjust = 1, color = dark_gray, size = 3.0, lineheight = 1.25) +

  annotate("text", x = 2.95, y = 4.55, label = "Gene-level standardization",
           hjust = 0, color = dark_gray, fontface = "bold", size = 3.15) +
  annotate("text", x = 2.95, y = 4.08,
           label = "z(g,i) = [x(g,i) - gene mean] / gene SD",
           hjust = 0, color = navy, size = 3.05) +
  annotate("text", x = 2.95, y = 3.35,
           label = "Module(k,i) = mean available gene z-scores",
           hjust = 0, color = dark_gray, size = 2.95) +
  annotate("text", x = 2.95, y = 2.82,
           label = "Summary(i) = mean z(HLA-II,\nCXCL9/FOLR2, DC/HLA)",
           hjust = 0, vjust = 0.5, color = dark_gray, size = 2.9,
           lineheight = 1.12) +
  annotate("segment", x = 2.95, xend = 5.05, y = 2.2, yend = 2.2,
           color = "#CDD3D7", linewidth = 0.5) +
  annotate("text", x = 2.95, y = 1.78,
           label = "Shared HLA genes disclosed",
           hjust = 0, color = orange, fontface = "bold", size = 2.9) +
  annotate("text", x = 2.95, y = 1.31,
           label = "Coverage checked; weights fixed\nafter refinement",
           hjust = 0, vjust = 0.5, color = dark_gray, size = 2.75,
           lineheight = 1.12) +

  annotate("point", x = 5.72, y = c(4.55, 3.78, 3.01, 2.24, 1.47),
           size = 3.2, shape = 21, stroke = 0.7,
           fill = c(blue, navy, teal, teal, orange), color = "white") +
  annotate("text", x = 5.95, y = c(4.55, 3.78, 3.01, 2.24, 1.47),
           label = c("GSE194040 / I-SPY2", "GSE25066 RMA",
                     "Five NAC cohorts + REML", "METABRIC", "GSE266919"),
           hjust = 0, vjust = -0.25, color = dark_gray, fontface = "bold", size = 2.9) +
  annotate("text", x = 5.95, y = c(4.55, 3.78, 3.01, 2.24, 1.47),
           label = c("internal association", "independent chemotherapy test",
                     "external synthesis", "recurrence biology", "cell-source localization"),
           hjust = 0, vjust = 1.65, color = c(blue, navy, teal, teal, orange), size = 2.55) +

  annotate("text", x = 9.12, y = 4.57, label = "Supported interpretation",
           hjust = 0, color = teal, fontface = "bold", size = 3.2) +
  annotate("text", x = 9.12, y = 4.14,
           label = "Context-dependent antigen-presentation state\nChemotherapy-context heterogeneity\nImmune-compartment localization",
           hjust = 0, vjust = 1, color = dark_gray, size = 2.9, lineheight = 1.24) +
  annotate("text", x = 9.12, y = 2.45, label = "Not established",
           hjust = 0, color = orange, fontface = "bold", size = 3.2) +
  annotate("text", x = 9.12, y = 2.02,
           label = "Universal pCR prediction\nIndependent component dimensions\nClinical utility",
           hjust = 0, vjust = 1, color = dark_gray, size = 2.9, lineheight = 1.24) +

  annotate("text", x = 0.32, y = 0.48, label = "Analysis chronology",
           hjust = 0, color = navy, fontface = "bold", size = 2.95) +
  annotate("text", x = 2.0, y = 0.48,
           label = "internal association  |  independent negative context  |  external synthesis  |  biological localization",
           hjust = 0, color = dark_gray, size = 2.8) +
  coord_cartesian(xlim = c(0, 12), ylim = c(0.15, 6.0), clip = "off") +
  theme_void(base_family = "Helvetica") +
  theme(plot.margin = margin(10, 14, 10, 14))

save_all(p1, "Figure1_experimental_technical_design_v2_6", 12.0, 6.2)

# Figure 2: response, proliferation, and recurrence forests
f2 <- read_tsv(file.path(src, "Figure2_effects_source_data.tsv"), show_col_types = FALSE) %>%
  mutate(
    role = factor(role, levels = c("internal association", "independent chemotherapy context",
                                  "external cohort", "external meta-analysis", "sensitivity")),
    effect_text = sprintf("%.2f (%.2f-%.2f)", estimate, conf_low, conf_high)
  )

role_fill <- c(
  "internal association" = blue,
  "independent chemotherapy context" = navy,
  "external cohort" = teal,
  "external meta-analysis" = dark_gray,
  "sensitivity" = "white"
)
role_color <- c(
  "internal association" = blue,
  "independent chemotherapy context" = navy,
  "external cohort" = teal,
  "external meta-analysis" = dark_gray,
  "sensitivity" = navy
)

forest_panel <- function(dat, title, xlab, limits, breaks, label_x) {
  dat <- dat %>% mutate(display_label = factor(display_label, levels = rev(unique(display_label))))
  ggplot(dat, aes(x = estimate, y = display_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = mid_gray, linewidth = 0.5) +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high, color = role),
                  orientation = "y", width = 0.16, linewidth = 0.65) +
    geom_point(aes(fill = role, color = role, shape = point_shape), size = 2.8, stroke = 0.75) +
    geom_text(aes(x = label_x, label = effect_text), hjust = 0, size = 2.45,
              color = dark_gray, family = "Helvetica") +
    scale_x_log10(limits = limits, breaks = breaks, labels = label_number(accuracy = 0.1)) +
    scale_fill_manual(values = role_fill, drop = FALSE) +
    scale_color_manual(values = role_color, drop = FALSE) +
    scale_shape_manual(values = c("circle" = 21, "diamond" = 23, "open" = 21)) +
    labs(title = title, x = xlab, y = NULL, fill = NULL, color = NULL, shape = NULL) +
    theme_pub(8.5) +
    theme(legend.position = "none", axis.text.y = element_text(size = 7.4))
}

p2a <- forest_panel(filter(f2, panel == "A_pCR_summary"),
                    "A  Antigen-presentation summary and pCR", "Odds ratio per 1-SD increase",
                    c(0.28, 5.6), c(0.3, 0.5, 1, 2, 4), 2.95)
p2b <- forest_panel(filter(f2, panel == "B_pCR_proliferation"),
                    "B  Proliferation comparator", "Odds ratio per 1-SD increase",
                    c(0.45, 7.2), c(0.5, 1, 2, 4), 4.35)
p2c <- forest_panel(filter(f2, panel == "C_recurrence"),
                    "C  Recurrence endpoints", "Hazard ratio per 1-SD increase",
                    c(0.45, 1.95), c(0.5, 0.75, 1, 1.5), 1.15)

p2 <- (p2a | (p2b / p2c)) +
  plot_layout(widths = c(1.45, 1)) +
  plot_annotation(
    caption = paste(
      "Blue: internal or independent chemotherapy context; teal: additional external cohorts;",
      "gray diamond: pooled estimate; hollow point: age- and stage-adjusted complete-case sensitivity."
    ),
    theme = theme(
      plot.caption = element_text(size = 8, color = dark_gray, hjust = 0),
      plot.margin = margin(4, 4, 4, 4)
    )
  )
save_all(p2, "Figure2_context_specific_effects_v2_6", 12.2, 8.2)

# Figure 3: score composition and correlation structure
overlap <- read_tsv(file.path(src, "Figure3_gene_overlap.tsv"), show_col_types = FALSE) %>%
  arrange(desc(n_genes)) %>%
  mutate(combination = factor(combination, levels = rev(combination)))

p3a <- ggplot(overlap, aes(n_genes, combination)) +
  geom_col(fill = navy, width = 0.66) +
  geom_text(aes(label = paste0(n_genes, ": ", genes)), hjust = -0.03, size = 2.35, color = dark_gray) +
  coord_cartesian(xlim = c(0, max(overlap$n_genes) * 2.45), clip = "off") +
  labs(title = "A  Exact gene-set membership patterns", x = "Number of genes", y = NULL) +
  theme_pub(8.5) + theme(axis.text.y = element_text(size = 7.2))

corr <- read_tsv(file.path(src, "Figure3_component_correlations.tsv"), show_col_types = FALSE)
p3b <- ggplot(corr, aes(score_2_label, score_1_label, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), fontface = "bold", size = 3.1,
            color = "white") +
  scale_fill_gradient2(low = "#B3584A", mid = "white", high = navy, midpoint = 0,
                       limits = c(-1, 1), name = "Spearman\nrho") +
  labs(title = "B  Component correlation in GSE194040", x = NULL, y = NULL) +
  theme_pub(8.5) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1), panel.grid = element_blank())

cross <- read_tsv(file.path(src, "Figure3_cross_cohort_correlations.tsv"), show_col_types = FALSE)
p3c <- ggplot(cross, aes(score_label, dataset_label, fill = spearman)) +
  geom_tile(color = "white", linewidth = 0.65) +
  geom_text(aes(label = sprintf("%.2f", spearman)), size = 2.55,
            color = ifelse(abs(cross$spearman) > 0.62, "white", dark_gray)) +
  scale_fill_gradient2(low = "#B3584A", mid = "white", high = navy, midpoint = 0,
                       limits = c(-1, 1), name = "Spearman\nrho") +
  labs(title = "C  Correlation with contextual scores", x = NULL, y = NULL) +
  theme_pub(8.5) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1), panel.grid = element_blank())

adj <- read_tsv(file.path(src, "Figure3_proliferation_adjustment.tsv"), show_col_types = FALSE) %>%
  mutate(dataset = factor(dataset, levels = rev(unique(dataset))),
         model = factor(model, levels = c("Summary only", "Summary + proliferation")))
p3d <- ggplot(adj, aes(estimate, dataset, group = dataset)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = mid_gray) +
  geom_line(aes(group = dataset), color = "#B8BEC2", linewidth = 0.55) +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high, color = model),
                orientation = "y", width = 0.12, linewidth = 0.55) +
  geom_point(aes(fill = model, color = model), shape = 21, size = 2.5, stroke = 0.65) +
  scale_x_log10(limits = c(0.3, 3.2), breaks = c(0.5, 1, 2, 3)) +
  scale_color_manual(values = c("Summary only" = navy, "Summary + proliferation" = orange)) +
  scale_fill_manual(values = c("Summary only" = navy, "Summary + proliferation" = "white")) +
  labs(title = "D  Effect after proliferation adjustment", x = "Odds ratio per 1-SD increase", y = NULL,
       color = NULL, fill = NULL) +
  theme_pub(8.5) + theme(legend.position = "bottom", axis.text.y = element_text(size = 7.2))

p3 <- (p3a | p3b) / (p3c | p3d) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
save_all(p3, "Figure3_summary_score_transparency_v2_5", 11.8, 8.2)

# Figure 4: actual single-cell subcluster scores
f4 <- read_tsv(file.path(src, "Figure4_scRNA_source_data.tsv"), show_col_types = FALSE)
main4 <- f4 %>%
  filter(source_interpretation == "Selected coherent APC/myeloid/B-cell subcluster",
         !grepl("MAST", subcluster),
         source_compartment %in% c("Myeloid", "Bcell")) %>%
  mutate(
    module_label = factor(module_label,
                          levels = c("HLA-II", "CXCL9/FOLR2 macrophage", "Dendritic-cell/HLA")),
    source_compartment = recode(source_compartment, Bcell = "B cell"),
    subcluster = factor(subcluster, levels = rev(unique(subcluster[order(source_compartment, -score_z_within_module)])))
  )

p4 <- ggplot(main4, aes(module_label, subcluster)) +
  geom_point(aes(size = n_patients, color = score_z_within_module), alpha = 0.92) +
  facet_grid(source_compartment ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_color_gradient2(low = "#C76752", mid = "#F4F4F2", high = navy, midpoint = 0,
                        name = "Within-module\nstandardized score") +
  scale_size_continuous(range = c(2.5, 7), breaks = c(30, 35, 40, 43), name = "Patients") +
  labs(x = NULL, y = NULL) +
  theme_pub(9) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"),
    strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold", color = navy),
    legend.position = "right"
  )
save_all(p4, "Figure4_GSE266919_cell_source_v2_5", 8.7, 7.5)

supp4 <- f4 %>%
  mutate(
    module_label = factor(module_label,
                          levels = c("HLA-II", "CXCL9/FOLR2 macrophage", "Dendritic-cell/HLA")),
    source_compartment = recode(source_compartment, Bcell = "B cell", CD4Tcell = "CD4 T cell", CD8Tcell = "CD8 T cell"),
    subcluster = factor(subcluster, levels = rev(unique(subcluster[order(source_compartment, -score_z_within_module)])))
  )

ps <- ggplot(supp4, aes(module_label, subcluster)) +
  geom_point(aes(size = n_patients, color = score_z_within_module,
                 alpha = source_interpretation != "Context only")) +
  facet_grid(source_compartment ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_color_gradient2(low = "#C76752", mid = "#F4F4F2", high = navy, midpoint = 0,
                        name = "Within-module\nstandardized score") +
  scale_size_continuous(range = c(2.2, 6), name = "Patients") +
  scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.38), guide = "none") +
  labs(x = NULL, y = NULL) + theme_pub(8.2) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold"))
save_all(ps, "Supplementary_Figure1_GSE266919_detailed_source_v2_5", 8.9, 9.8, supplementary = TRUE)

cat("Rendered Figure 1-4 and Supplementary Figure 1\n")
