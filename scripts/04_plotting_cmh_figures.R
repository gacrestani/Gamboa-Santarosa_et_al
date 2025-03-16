# plotting_cmh_figures.R
# Plots all the cmh figures possible. Later I can choose which ones are the best.
# Note: cmh plotting functions work well for when you're plotting the same sample in two different time points
# If you want to plot one sample vs other (i.e. OBO_20 vs OB_20), it may not work!
#
# inputs: snp_table.rds (x4), cmh_pvals.rds, perm_pvals.csv
# outputs: manhattan plots, in accordance with the following format:
#
# ./results/figures/cmh
# │   │   ├── cmh_adapted_B.png
# │   │   ├── cmh_adapted_nBO.png
# │   │   ├── cmh_adapted_nB.png
# │   │   ├── cmh_adapted_OBO_OB_O_piled.png
# │   │   ├── cmh_adapted_OBO.png
# │   │   ├── cmh_adapted_OB.png
# │   │   └── cmh_adapted_O.png
# │       ├── cmh_classic_B.png
# │       ├── cmh_classic_nBO.png
# │       ├── cmh_classic_nB.png
# │       ├── cmh_classic_OBO_OB_O_piled.png
# │       ├── cmh_classic_OBO.png
# │       ├── cmh_classic_OB.png
# │       └── cmh_classic_O.png
# │   │   ├── cmh_adapted_fdr_B.png
# │   │   ├── cmh_adapted_fdr_nBO.png
# │   │   ├── cmh_adapted_fdr_nB.png
# │   │   ├── cmh_adapted_fdr_OBO_OB_O_piled.png
# │   │   ├── cmh_adapted_fdr_OBO.png
# │   │   ├── cmh_adapted_fdr_OB.png
# │   │   └── cmh_adapted_fdr_O.png
# │       ├── cmh_classic_fdr_B.png
# │       ├── cmh_classic_fdr_nBO.png
# │       ├── cmh_classic_fdr_nB.png
# │       ├── cmh_classic_fdr_OBO_OB_O_piled.png
# │       ├── cmh_classic_fdr_OBO.png
# │       ├── cmh_classic_fdr_OB.png
# │       └── cmh_classic_fdr_O.png
# │   │   ├── cmh_adapted_fdr_scaled_B.png
# │   │   ├── cmh_adapted_fdr_scaled_nBO.png
# │   │   ├── cmh_adapted_fdr_scaled_nB.png
# │   │   ├── cmh_adapted_fdr_scaled_OBO_OB_O_piled.png
# │   │   ├── cmh_adapted_fdr_scaled_OBO.png
# │   │   ├── cmh_adapted_fdr_scaled_OB.png
# │   │   └── cmh_adapted_fdr_scaled_O.png
# │       ├── cmh_classic_fdr_scaled_B.png
# │       ├── cmh_classic_fdr_scaled_nBO.png
# │       ├── cmh_classic_fdr_scaled_nB.png
# │       ├── cmh_classic_fdr_scaled_OBO_OB_O_piled.png
# │       ├── cmh_classic_fdr_scaled_OBO.png
# │       ├── cmh_classic_fdr_scaled_OB.png
# │       └── cmh_classic_fdr_scaled_O.png
#     │   ├── cmh_adapted_scaled_B.png
#     │   ├── cmh_adapted_scaled_nBO.png
#     │   ├── cmh_adapted_scaled_nB.png
#     │   ├── cmh_adapted_scaled_OBO_OB_O_piled.png
#     │   ├── cmh_adapted_scaled_OBO.png
#     │   ├── cmh_adapted_scaled_OB.png
#     │   └── cmh_adapted_scaled_O.png
#         ├── cmh_classic_scaled_B.png
#         ├── cmh_classic_scaled_nBO.png
#         ├── cmh_classic_scaled_nB.png
#         ├── cmh_classic_scaled_OBO_OB_O_piled.png
#         ├── cmh_classic_scaled_OBO.png
#         ├── cmh_classic_scaled_OB.png
#         └── cmh_classic_scaled_O.png

source("scripts/functions.R")

# 0 Initializing ---------------------------------------------------------------
# Read snp tables
snp_table_shahrestani <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

snp_table_regimes <-
  readRDS("data/processed/processed_snps_abcd_regimes.rds")

snp_table_shahrestani_scaled <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani_scaled.rds")

snp_table_regimes_scaled <-
  readRDS("data/processed/processed_snps_abcd_regimes_scaled.rds")

cmh_pvals <- readRDS("results/cmh_pvals.rds")
perm_pvals <- fread("results/perm_pvals.csv")

# Parameters
y_limit_up <- 220 # manhattan plot y-axis upper limit
width      <- 7740/2 # manhattan plot width
height     <- 1440/2 # manhattan plot height

treatments <- c("OBO", "OB", "nBO", "nB", "O", "B")
gen2 <- c("20", "20", "56", "56", "20", "56")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 1 Classic CMH ----------------------------------------------------------------
for (treatment in treatments) {
  plot <-
    GetManhattanPlot(
      my_dataframe = cmh_pvals,
      Y = -log10(cmh_pvals[[paste0("cmh_classic_", treatment)]]),
      permutation_pvals = perm_pvals[[treatments[i]]],
      percentage_significance = FALSE,
      title = paste0("Classical CMH test: ", treatments[i], " gen01 vs ", treatments[i], " gen", gen2[i]),
      x_label = TRUE,
      y_label = "-log10(p-value)",
      palette = "blue",
      y_limit_up = y_limit_up,
      y_limit_down = 0)
  
  ggsave(
    paste0("results/figures/cmh/cmh_classic_", treatments[i], ".png"),
    plot = plot,
    width = width,
    height = height,
    bg = "white",
    units = "px")
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
# 2 Adapted CMH ----------------------------------------------------------------
for (treatment in treatments) {
  plot <-
    GetManhattanPlot(
      my_dataframe = cmh_pvals,
      Y = -log10(cmh_pvals[[paste0("cmh_adapted_", treatments[i])]]),
      permutation_pvals = NULL,
      percentage_significance = TRUE,
      title = paste0("Adapted CMH test: ", treatments[i], " gen01 vs ", treatments[i], " gen", gen2[i]),
      x_label = TRUE,
      y_label = "-log10(p-value)",
      palette = "blue",
      y_limit_up = y_limit_up,
      y_limit_down = 0)
  
  ggsave(
    paste0("results/figures/cmh/cmh_adapted_", treatments[i], ".png"),
    plot = plot,
    width = width,
    height = height,
    bg = "white",
    units = "px")
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 3 Classic CMH Scaled ---------------------------------------------------------
for (treatment in treatments) {
  plot <-
    GetManhattanPlot(
      my_dataframe = cmh_pvals,
      Y = -log10(cmh_pvals[[paste0("cmh_classic_scaled_", treatments[i])]]),
      permutation_pvals = NULL,
      percentage_significance = TRUE,
      title = paste0("Classical CMH test: ", treatments[i], " gen01 vs ", treatments[i], " gen", gen2[i], ", scaled coverage."),
      x_label = TRUE,
      y_label = "-log10(p-value)",
      palette = "blue",
      y_limit_up = y_limit_up,
      y_limit_down = 0)
  
  ggsave(
    paste0("results/figures/cmh/cmh_classic_scaled_", treatments[i], ".png"),
    plot = plot,
    width = width,
    height = height,
    bg = "white",
    units = "px")
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 3 Adapted CMH Scaled ---------------------------------------------------------
for (treatment in treatments) {
  plot <-
    GetManhattanPlot(
      my_dataframe = cmh_pvals,
      Y = -log10(cmh_pvals[[paste0("cmh_adapted_scaled_", treatments[i])]]),
      permutation_pvals = NULL,
      percentage_significance = TRUE,
      title = paste0("Adapted CMH test: ", treatments[i], " gen01 vs ", treatments[i], " gen", gen2[i], ", scaled coverage."),
      x_label = TRUE,
      y_label = "-log10(p-value)",
      palette = "blue",
      y_limit_up = y_limit_up,
      y_limit_down = 0)
  
  ggsave(
    paste0("results/figures/cmh/cmh_adapted_scaled_", treatments[i], ".png"),
    plot = plot,
    width = width,
    height = height,
    bg = "white",
    units = "px")
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 4 Classic CMH FDR ------------------------------------------------------------
for (treatment in treatments) {
  plot <-
    GetManhattanPlot(
      my_dataframe = cmh_pvals,
      Y = -log10(p.adjust(cmh_pvals[[paste0("cmh_classic_scaled_", treatments[i])]], method = "BH")),
      permutation_pvals = NULL,
      percentage_significance = TRUE,
      title = paste0("Classic CMH test: ", treatments[i], " gen01 vs ", treatments[i], " gen", gen2[i], "FDR corrected."),
      x_label = TRUE,
      y_label = "-log10(p-value)",
      palette = "blue",
      y_limit_up = y_limit_up,
      y_limit_down = 0)
  
  ggsave(
    paste0("results/figures/cmh/cmh_classic_fdr_", treatments[i], ".png"),
    plot = plot,
    width = width,
    height = height,
    bg = "white",
    units = "px")
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 5 Adapted CMH FDR ------------------------------------------------------------
for (treatment in treatments) {
  plot <-
    GetManhattanPlot(
      my_dataframe = cmh_pvals,
      Y = -log10(p.adjust(cmh_pvals[[paste0("cmh_classic_scaled_", treatments[i])]], method = "BH")),
      permutation_pvals = NULL,
      percentage_significance = TRUE,
      title = paste0("Adapted CMH test: ", treatments[i], " gen01 vs ", treatments[i], " gen", gen2[i], "FDR corrected."),
      x_label = TRUE,
      y_label = "-log10(p-value)",
      palette = "blue",
      y_limit_up = y_limit_up,
      y_limit_down = 0)
  
  ggsave(
    paste0("results/figures/cmh/cmh_adapted_fdr_", treatments[i], ".png"),
    plot = plot,
    width = width,
    height = height,
    bg = "white",
    units = "px")
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# 6 Classic CMH Scaled FDR -----------------------------------------------------
for (treatment in treatments) {
  plot <-
    GetManhattanPlot(
      my_dataframe = cmh_pvals,
      Y = -log10(p.adjust(cmh_pvals[[paste0("cmh_classic_scaled_", treatments[i])]], method = "BH")),
      permutation_pvals = NULL,
      percentage_significance = TRUE,
      title = paste0("Classic CMH test: ", treatments[i], " gen01 vs ", treatments[i], " gen", gen2[i], ", scaled coverage, FDR corrected."),
      x_label = TRUE,
      y_label = "-log10(p-value)",
      palette = "blue",
      y_limit_up = y_limit_up,
      y_limit_down = 0)
  
  ggsave(
    paste0("results/figures/cmh/cmh_classic_scaled_fdr_", treatments[i], ".png"),
    plot = plot,
    width = width,
    height = height,
    bg = "white",
    units = "px")
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 7 Adapted CMH Scaled FDR -----------------------------------------------------
for (treatment in treatments) {
  plot <-
    GetManhattanPlot(
      my_dataframe = cmh_pvals,
      Y = -log10(p.adjust(cmh_pvals[[paste0("cmh_classic_scaled_", treatments[i])]], method = "BH")),
      permutation_pvals = NULL,
      percentage_significance = TRUE,
      title = paste0("Adapted CMH test: ", treatments[i], " gen01 vs ", treatments[i], " gen", gen2[i], ", scaled coverage, FDR corrected."),
      x_label = TRUE,
      y_label = "-log10(p-value)",
      palette = "blue",
      y_limit_up = y_limit_up,
      y_limit_down = 0)
  
  ggsave(
    paste0("results/figures/cmh/cmh_adapted_scaled_fdr_", treatments[i], ".png"),
    plot = plot,
    width = width,
    height = height,
    bg = "white",
    units = "px")
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 8 Grid plots -----------------------------------------------------------------

layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)

grid_plot_cmh_classic_OBO <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_obo01_vs_obo20),
    permutation_pvals = perm_pvals$obo,
    percentage_significance = FALSE,
    title = "Classical CMH test: OBO gen01 vs OBO gen20",
    x_label = FALSE,
    y_label = NULL,
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

grid_plot_cmh_classic_OB <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_ob01_vs_ob20),
    permutation_pvals = perm_pvals$ob,
    percentage_significance = FALSE,
    title = "Classical CMH test: OB gen01 vs OB gen20",
    x_label = FALSE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

grid_plot_cmh_classic_O <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_o01_vs_o20),
    permutation_pvals = perm_pvals$o,
    percentage_significance = FALSE,
    title = "Classical CMH test: O gen01 vs O gen20",
    x_label = TRUE,
    y_label = NULL,
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

png(
  filename = "results/figures/cmh_crude/classic/cmh_classic_OBO_OB_O_piled.png",
  width = 1800,
  height = 900)

grid.arrange(
  grid_plot_cmh_classic_OBO,
  grid_plot_cmh_classic_OB,
  grid_plot_cmh_classic_O,
  layout_matrix = layout)

dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~