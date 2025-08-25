# plotting_cmh_figures.R
# Plots all the cmh figures possible. Later I can choose which ones are the best.
# Note: cmh plotting functions work well for when you're plotting the same sample in two different time points
# If you want to plot one sample vs other (i.e. OBO_20 vs OB_20), it may not work!
#
# inputs: snp_table.rds (x4), cmh_pvals.rds, perm_pvals.csv
# outputs: manhattan plots

source("scripts/functions.R")

# 0 Initializing ---------------------------------------------------------------
# Read snp tables
# snp_table_shahrestani <- 
#   readRDS("data/processed/processed_snps_abcd_shahrestani.rds")
# 
# snp_table_regimes <-
#   readRDS("data/processed/processed_snps_abcd_regimes.rds")

snp_table_shahrestani_scaled <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani_scaled.rds")

snp_table_regimes_scaled <-
  readRDS("data/processed/processed_snps_abcd_regimes_scaled.rds")

cmh_pvals <- readRDS("results/cmh_pvals.rds")
perm_pvals <- fread("results/perm_pvals.csv")

# Parameters
y_limit_up <- 220 # manhattan plot y-axis upper limit

parameters <- data.frame(matrix(ncol = 4, nrow = 48))
colnames(parameters) <-
  c(
    "treatments",
    "gen2", 
    "pvals",
    "title"
    )

treatments <- c("OBO", "OB", "nBO", "nB", "O",  "B")
parameters$treatments <- rep(treatments, times = 4)

gen2 <-       c("20",  "20", "56",  "56", "20", "56")
parameters$gen2 <- rep(gen2, times = 4)

parameters$pvals <- colnames(cmh_pvals[-c(1,2)])

titles <-
  c(
    "Classic CMH test:",
    "Adapted CMH test:", 
    "Classic CMH test, scaled:", 
    "Adapted CMH test, scaled:", 
    "Classic CMH test, FDR corrected:",
    "Adapted CMH test, FDR corrected:", 
    "Classic CMH test, scaled, FDR corrected:", 
    "Adapted CMH test, scaled, FDR corrected:" 
    )

parameters$title <- rep(titles, each = 6)
parameters$title <- paste(parameters$title, parameters$treatments, "gen 01 vs", parameters$treatments, "gen", parameters$gen2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mapply(function(pvals, title, filename) {
  plot <- GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals[[pvals]]),
    percentage_significance = TRUE,
    title = title,
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0
  )
  
  ggsave(
    filename = paste0("results/figures/cmh/", pvals, ".png"),
    plot = plot,
    width = width,
    height = height,
    bg = "white",
    units = "px"
  )
  
  print("Concluded", parameters$title)
}, parameters$pvals, parameters$title, SIMPLIFY = FALSE)

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


obo <- GetManhattanPlot(
  my_dataframe = cmh_pvals,
  Y = -log10(cmh_pvals$cmh_adapted_fdr_obo01_vs_obo20_scaled),
  #permutation_pvals = perm_pvals$o,
  percentage_significance = TRUE,
  title = "OBO",
  x_label = FALSE,
  y_label = "",
  palette = "blue",
  y_limit_up = y_limit_up,
  y_limit_down = 0)

ob <- GetManhattanPlot(
  my_dataframe = cmh_pvals,
  Y = -log10(cmh_pvals$cmh_adapted_fdr_ob01_vs_ob20_scaled),
  #permutation_pvals = perm_pvals$o,
  percentage_significance = TRUE,
  title = "OB",
  x_label = FALSE,
  y_label = "-log10(q-values)",
  palette = "blue",
  y_limit_up = y_limit_up,
  y_limit_down = 0)

o <- GetManhattanPlot(
  my_dataframe = cmh_pvals,
  Y = -log10(cmh_pvals$cmh_adapted_o01_vs_o20),
  #permutation_pvals = perm_pvals$o,
  percentage_significance = TRUE,
  title = "O",
  x_label = TRUE,
  y_label = "",
  palette = "blue",
  y_limit_up = y_limit_up,
  y_limit_down = 0)


# ggsave("results/figures/cmh/cmh_fdr_scaled/adapted/obo.png", obo, width = 8, height = 3, units = "in", dpi = 450)
# ggsave("results/figures/cmh/cmh_fdr_scaled/adapted/ob.png", ob, width = 8, height = 3, units = "in", dpi = 450)
# ggsave("results/figures/cmh/cmh_fdr_scaled/adapted/o.png", o, width = 8, height = 3, units = "in", dpi = 450)

grid_plot <- grid.arrange(obo,
                          ob,
                          o,
                          nrow = 3)

ggsave("results/figures/cmh/cmh_fdr_scaled/adapted/obo_ob_o.png", grid_plot, width = 7.5, height = 9, units = "in", dpi = 600)

nbo <- GetManhattanPlot(
  my_dataframe = cmh_pvals,
  Y = -log10(cmh_pvals$cmh_adapted_fdr_nbo01_vs_nbo56_scaled),
  percentage_significance = TRUE,
  title = "nBO",
  x_label = FALSE,
  y_label = "",
  palette = "red",
  y_limit_up = y_limit_up,
  y_limit_down = 0
) 

nb <- GetManhattanPlot(
  my_dataframe = cmh_pvals,
  Y = -log10(cmh_pvals$cmh_adapted_fdr_nb01_vs_nb56_scaled),
  percentage_significance = TRUE,
  title = "nB",
  x_label = FALSE,
  y_label = "-log10(q-value)",
  palette = "red",
  y_limit_up = y_limit_up,
  y_limit_down = 0
) 

b <- GetManhattanPlot(
  my_dataframe = cmh_pvals,
  Y = -log10(cmh_pvals$cmh_adapted_fdr_b01_vs_b56_scaled),
  percentage_significance = TRUE,
  title = "B",
  x_label = TRUE,
  y_label = "",
  palette = "red",
  y_limit_up = y_limit_up,
  y_limit_down = 0
) 

# ggsave("results/figures/cmh/cmh_fdr_scaled/adapted/b.png", b, width = 8, height = 3, units = "in", dpi = 450)
# ggsave("results/figures/cmh/cmh_fdr_scaled/adapted/b.png", b, width = 8, height = 3, units = "in", dpi = 450)
# ggsave("results/figures/cmh/cmh_fdr_scaled/adapted/b.png", b, width = 8, height = 3, units = "in", dpi = 450)

grid_plot <- grid.arrange(nbo,
                          nb,
                          b,
                          nrow = 3)

ggsave("results/figures/cmh/cmh_fdr_scaled/adapted/nbo_nb_b.png", grid_plot, width = 7.5, height = 9, units = "in", dpi = 600)










# B56 vs O20

p_vals_ob <- ClassicalCmhTest(
  snp_table_shahrestani_scaled,
  treatment1 = "OBO",
  gen1 = "20",
  treatment2 = "nBO",
  gen2 = "56"
)

p_vals_b <- ClassicalCmhTest(
  snp_table_shahrestani_scaled,
  treatment1 = "OB",
  gen1 = "20",
  treatment2 = "nB",
  gen2 = "56"
)

p_vals_b_ob <- ClassicalCmhTest(
  snp_table_regimes_scaled,
  treatment1 = "O",
  gen1 = "20",
  treatment2 = "B",
  gen2 = "56"
)

p_vals_adj_ob <- p.adjust(p_vals_ob, method = "fdr")
p_vals_adj_b <- p.adjust(p_vals_b, method = "fdr")
p_vals_adj_b_ob <- p.adjust(p_vals_b_ob, method = "fdr")

ob <- GetManhattanPlot(
  my_dataframe = cmh_pvals,
  Y = -log10(p_vals_adj_ob),
  percentage_significance = TRUE,
  title = "OBO20 vs nBO55",
  x_label = FALSE,
  y_label = "",
  palette = "grey",
  y_limit_up = y_limit_up,
  y_limit_down = 0
) 


b <- GetManhattanPlot(
  my_dataframe = cmh_pvals,
  Y = -log10(p_vals_adj_b),
  percentage_significance = TRUE,
  title = "OB20 vs nB55",
  x_label = FALSE,
  y_label = "-log10(q-value)",
  palette = "grey",
  y_limit_up = y_limit_up,
  y_limit_down = 0
) 

b_ob <- GetManhattanPlot(
  my_dataframe = cmh_pvals,
  Y = -log10(p_vals_adj_b_ob),
  percentage_significance = TRUE,
  title = "O20 vs B55",
  x_label = TRUE,
  y_label = "",
  palette = "grey",
  y_limit_up = y_limit_up,
  y_limit_down = 0
) 

grid_plot <- grid.arrange(ob,
                          b,
                          b_ob,
                          nrow = 3)

ggsave("results/figures/cmh/cmh_fdr_scaled/adapted/o20_vs_b56_classic.png", grid_plot, width = 7.5, height = 9, units = "in", dpi = 600)

