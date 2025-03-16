
# 4 - Plotting the results =====================================================
cmh_pvals <- readRDS("results/cmh_pvals.rds")
perm_pvals <- fread("results/perm_pvals.csv")



# 4.1 - Crude ==================================================================
# 4.1.1 - Classic CMH ==========================================================
plot_cmh_classic_OBO <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_obo01_vs_obo20),
    permutation_pvals = perm_pvals$obo,
    title = "Classical CMH test: OBO gen01 vs OBO gen20",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/classic/cmh_classic_OBO.png",
  plot = plot_cmh_classic_OBO,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_classic_OB <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_ob01_vs_ob20),
    permutation_pvals = perm_pvals$ob,
    title = "Classical CMH test: OB gen01 vs OB gen20",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/classic/cmh_classic_OB.png",
  plot = plot_cmh_classic_OB,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_classic_nBO <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_nbo01_vs_nbo56),
    permutation_pvals = perm_pvals$nbo,
    title = "Classical CMH test: nBO gen01 vs nBO gen56",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "red",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/classic/cmh_classic_nBO.png",
  plot = plot_cmh_classic_nBO,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_classic_nB <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_nb01_vs_nb56),
    permutation_pvals = perm_pvals$nb,
    title = "Classical CMH test: nB gen01 vs nB gen56",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "red",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/classic/cmh_classic_nB.png",
  plot = plot_cmh_classic_nB,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_classic_O <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_o01_vs_o20),
    permutation_pvals = perm_pvals$o,
    title = "Classical CMH test: O gen01 vs O gen20",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/classic/cmh_classic_O.png",
  plot = plot_cmh_classic_O,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_classic_B <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_classic_b01_vs_b56),
    permutation_pvals = perm_pvals$b,
    title = "Classical CMH test: B gen01 vs B gen56",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "red",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/classic/cmh_classic_B.png",
  plot = plot_cmh_classic_B,
  width = width,
  height = height,
  bg = "white",
  units = "px")

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

# 4.1.2 - Adapted CMH ==========================================================
plot_cmh_adapted_OBO <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_obo01_vs_obo20),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test: OBO gen01 vs OBO gen20",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/adapted/cmh_adapted_OBO.png",
  plot = plot_cmh_adapted_OBO,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_adapted_OB <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_ob01_vs_ob20),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test: OB gen01 vs OB gen20",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/adapted/cmh_adapted_OB.png",
  plot = plot_cmh_adapted_OB,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_adapted_nBO <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_nbo01_vs_nbo56),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test: nBO gen01 vs nBO gen56",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "red",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/adapted/cmh_adapted_nBO.png",
  plot = plot_cmh_adapted_nBO,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_adapted_nB <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_nb01_vs_nb56),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test: nB gen01 vs nB gen56",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "red",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/adapted/cmh_adapted_nB.png",
  plot = plot_cmh_adapted_nB,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_adapted_O <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_o01_vs_o20),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test: O gen01 vs O gen20",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/adapted/cmh_adapted_O.png",
  plot = plot_cmh_adapted_O,
  width = width,
  height = height,
  bg = "white",
  units = "px")

plot_cmh_adapted_B <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_b01_vs_b56),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test: B gen01 vs B gen56",
    x_label = TRUE,
    y_label = "-log10(p-value)",
    palette = "red",
    y_limit_up = y_limit_up,
    y_limit_down = 0)
ggsave(
  "results/figures/cmh_crude/adapted/cmh_adapted_B.png",
  plot = plot_cmh_adapted_B,
  width = width,
  height = height,
  bg = "white",
  units = "px")

layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)

grid_plot_cmh_adapted_OBO <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_obo01_vs_obo20),
    permutation_pvals = obo_adapted$V1,
    percentage_significance = TRUE,
    title = "Adapted CMH test: OBO gen01 vs OBO gen20",
    x_label = FALSE,
    y_label = NULL,
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

grid_plot_cmh_adapted_OB <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_ob01_vs_ob20),
    permutation_pvals = ob_adapted$V1,
    percentage_significance = TRUE,
    title = "Adapted CMH test: OB gen01 vs OB gen20",
    x_label = FALSE,
    y_label = "-log10(p-value)",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

grid_plot_cmh_adapted_O <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(cmh_pvals$cmh_adapted_o01_vs_o20),
    permutation_pvals = o_adapted$V1,
    percentage_significance = TRUE,
    title = "Adapted CMH test: O gen01 vs O gen20",
    x_label = TRUE,
    y_label = NULL,
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

png(
  filename = "results/figures/cmh_crude/adapted/cmh_adapted_OBO_OB_O_piled.png",
  width = 1800,
  height = 900)

grid.arrange(
  grid_plot_cmh_adapted_OBO,
  grid_plot_cmh_adapted_OB,
  grid_plot_cmh_adapted_O,
  layout_matrix = layout)

dev.off()

# 4.2 - Scaled data ============================================================
# 4.2.1 - Classic CMH Scaled ===================================================
plot_cmh_classic_OBO_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_obo01_vs_obo20_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: OBO gen01 vs OBO gen20",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/classic/cmh_classic_OBO_scaled.png",
       plot = plot_cmh_classic_OBO_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_OB_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_ob01_vs_ob20_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: OB gen01 vs OB gen20",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/classic/cmh_classic_OB_scaled.png",
       plot = plot_cmh_classic_OB_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_nBO_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_nbo01_vs_nbo56_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: nBO gen01 vs nBO gen56",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/classic/cmh_classic_nBO_scaled.png",
       plot = plot_cmh_classic_nBO_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_nB_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_nb01_vs_nb56_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: nB gen01 vs nB gen56",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/classic/cmh_classic_nB_scaled.png",
       plot = plot_cmh_classic_nB_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_O_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_o01_vs_o20_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: O gen01 vs O gen20",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/classic/cmh_classic_O_scaled.png",
       plot = plot_cmh_classic_O_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_B_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_b01_vs_b56_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: B gen01 vs B gen56",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/classic/cmh_classic_B_scaled.png",
       plot = plot_cmh_classic_B_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

# Grid arrange  
layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)

grid_plot_cmh_classic_OBO_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_obo01_vs_obo20_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: OBO gen01 vs OBO gen20",
           x_label = FALSE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

grid_plot_cmh_classic_OB_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_ob01_vs_ob20_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: OB gen01 vs OB gen20",
           x_label = FALSE,
           y_label = "-log10(p-value)",
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

grid_plot_cmh_classic_O_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_classic_o01_vs_o20_scaled),
           permutation_pvals = NULL,
           title = "Classical CMH test, scaled: O gen01 vs O gen20",
           x_label = TRUE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

# Save grid
png(filename = "results/figures/cmh_scaled/classic/cmh_classic_OBO_OB_O_scaled_piled.png",
    width = 1800,
    height = 900)

grid.arrange(grid_plot_cmh_classic_OBO_scaled,
             grid_plot_cmh_classic_OB_scaled,
             grid_plot_cmh_classic_O_scaled,
             layout_matrix = layout)

dev.off()

# 4.2.2 - Adapted CMH Scaled ===================================================
plot_cmh_adapted_OBO_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_obo01_vs_obo20_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: OBO gen01 vs OBO gen20",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/adapted/cmh_adapted_OBO_scaled.png",
       plot = plot_cmh_adapted_OBO_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_OB_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_ob01_vs_ob20_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: OB gen01 vs OB gen20",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/adapted/cmh_adapted_OB_scaled.png",
       plot = plot_cmh_adapted_OB_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_nBO_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_nbo01_vs_nbo56_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: nBO gen01 vs nBO gen56",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/adapted/cmh_adapted_nBO_scaled.png",
       plot = plot_cmh_adapted_nBO_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_nB_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_nb01_vs_nb56_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: nB gen01 vs nB gen56",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/adapted/cmh_adapted_nB_scaled.png",
       plot = plot_cmh_adapted_nB_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_O_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_o01_vs_o20_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: O gen01 vs O gen20",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/adapted/cmh_adapted_O_scaled.png",
       plot = plot_cmh_adapted_O_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_B_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_b01_vs_b56_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: B gen01 vs B gen56",
           x_label = TRUE,
           y_label = "-log10(p-value)",
           palette = "red",
           y_limit_up = y_limit_up,
           y_limit_down = 0)
ggsave("results/figures/cmh_scaled/adapted/cmh_adapted_B_scaled.png",
       plot = plot_cmh_adapted_B_scaled,
       width = width,
       height = height,
       bg = "white",
       units = "px")

# Grid arrange
layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)

grid_plot_cmh_adapted_OBO_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_obo01_vs_obo20_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: OBO gen01 vs OBO gen20",
           x_label = FALSE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

grid_plot_cmh_adapted_OB_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_ob01_vs_ob20_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: OB gen01 vs OB gen20",
           x_label = FALSE,
           y_label = "-log10(p-value)",
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

grid_plot_cmh_adapted_O_scaled <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
           Y = -log10(cmh_pvals$cmh_adapted_o01_vs_o20_scaled),
           permutation_pvals = NULL,
           percentage_significance = TRUE,
           title = "Adapted CMH test, scaled: O gen01 vs O gen20",
           x_label = TRUE,
           y_label = NULL,
           palette = "blue",
           y_limit_up = y_limit_up,
           y_limit_down = 0)

# Save grid
png(filename = "results/figures/cmh_scaled/adapted/cmh_adapted_OBO_OB_O_scaled_piled.png",
    width = 1800,
    height = 900)

grid.arrange(grid_plot_cmh_adapted_OBO_scaled,
             grid_plot_cmh_adapted_OB_scaled,
             grid_plot_cmh_adapted_O_scaled,
             layout_matrix = layout)

dev.off()

# 4.3 - FDR Corrected ==========================================================
# 4.3.1 - Classic CMH FDR ======================================================
plot_cmh_classic_OBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_obo01_vs_obo20, method = "BH")),
                   permutation_pvals = perm_pvals$obo,
                   title = "Classical CMH test, FDR corrected: OBO gen01 vs OBO gen20",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/classic/cmh_classic_OBO.png",
       plot = plot_cmh_classic_OBO,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_OB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_ob01_vs_ob20, method = "BH")),
                   permutation_pvals = perm_pvals$ob,
                   title = "Classical CMH test, FDR corrected: OB gen01 vs OB gen20",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/classic/cmh_classic_OB.png",
       plot = plot_cmh_classic_OB,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_nBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_nbo01_vs_nbo56, method = "BH")),
                   permutation_pvals = perm_pvals$nbo,
                   title = "Classical CMH test, FDR corrected: nBO gen01 vs nBO gen56",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "red",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/classic/cmh_classic_nBO.png",
       plot = plot_cmh_classic_nBO,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_nB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_nb01_vs_nb56, method = "BH")),
                   permutation_pvals = perm_pvals$nb,
                   title = "Classical CMH test, FDR corrected: nB gen01 vs nB gen56",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "red",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/classic/cmh_classic_nB.png",
       plot = plot_cmh_classic_nB,
       width = width,
       height = height,
       bg = "white",
       units = "px")


plot_cmh_classic_O <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_o01_vs_o20, method = "BH")),
                   permutation_pvals = perm_pvals$o,
                   title = "Classical CMH test, FDR corrected: O gen01 vs O gen20",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/classic/cmh_classic_O.png",
       plot = plot_cmh_classic_O,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_classic_B <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_b01_vs_b56, method = "BH")),
                   permutation_pvals = perm_pvals$b,
                   title = "Classical CMH test, FDR corrected: B gen01 vs B gen56",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "red",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/classic/cmh_classic_B.png",
       plot = plot_cmh_classic_B,
       width = width,
       height = height,
       bg = "white",
       units = "px")

# Grid arrange
layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)

grid_plot_cmh_classic_OBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_obo01_vs_obo20, method = "BH")),
                   permutation_pvals = perm_pvals$obo,
                   title = "Classical CMH test, FDR corrected: OBO gen01 vs OBO gen20",
                   x_label = FALSE,
                   y_label = NULL,
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

grid_plot_cmh_classic_OB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_ob01_vs_ob20, method = "BH")),
                   permutation_pvals = perm_pvals$ob,
                   title = "Classical CMH test, FDR corrected: OB gen01 vs OB gen20",
                   x_label = FALSE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

grid_plot_cmh_classic_O <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_classic_o01_vs_o20, method = "BH")),
                   permutation_pvals = NULL,
                   title = "Classical CMH test, FDR corrected: O gen01 vs O gen20",
                   x_label = TRUE,
                   y_label = NULL,
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

# Save grid
png(filename = "results/figures/cmh_fdr/classic/cmh_classic_OBO_OB_O_piled.png",
    width = 1800,
    height = 900)

grid.arrange(grid_plot_cmh_classic_OBO,
             grid_plot_cmh_classic_OB,
             grid_plot_cmh_classic_O,
             layout_matrix = layout)

dev.off()

# 4.3.2 - Adapted CMH FDR ======================================================
plot_cmh_adapted_OBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_obo01_vs_obo20, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: OBO gen01 vs OBO gen20",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/adapted/cmh_adapted_OBO.png",
       plot = plot_cmh_adapted_OBO,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_OB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_ob01_vs_ob20, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: OB gen01 vs OB gen20",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/adapted/cmh_adapted_OB.png",
       plot = plot_cmh_adapted_OB,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_nBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_nbo01_vs_nbo56, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: nBO gen01 vs nBO gen56",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "red",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/adapted/cmh_adapted_nBO.png",
       plot = plot_cmh_adapted_nBO,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_nB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_nb01_vs_nb56, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: nB gen01 vs nB gen56",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "red",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/adapted/cmh_adapted_nB.png",
       plot = plot_cmh_adapted_nB,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_O <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_o01_vs_o20, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: O gen01 vs O gen20",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/adapted/cmh_adapted_O.png",
       plot = plot_cmh_adapted_O,
       width = width,
       height = height,
       bg = "white",
       units = "px")

plot_cmh_adapted_B <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_b01_vs_b56, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: B gen01 vs B gen56",
                   x_label = TRUE,
                   y_label = "-log10(p-value)",
                   palette = "red",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)
ggsave("results/figures/cmh_fdr/adapted/cmh_adapted_B.png",
       plot = plot_cmh_adapted_B,
       width = width,
       height = height,
       bg = "white",
       units = "px")

# Grid arrange
layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)

grid_plot_cmh_adapted_OBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_obo01_vs_obo20, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = NULL,
                   x_label = FALSE,
                   y_label = NULL,
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

grid_plot_cmh_adapted_OB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_ob01_vs_ob20, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: OB gen01 vs OB gen20",
                   x_label = FALSE,
                   y_label = "-log10(p-value)",
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

grid_plot_cmh_adapted_O <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_o01_vs_o20, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = "Adapted CMH test, FDR corrected: O gen01 vs O gen20",
                   x_label = TRUE,
                   y_label = NULL,
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

# Save grid
png(filename = "results/figures/cmh_fdr/adapted/cmh_adapted_OBO_OB_O_piled.png",
    width = 1800,
    height = 900)

grid.arrange(grid_plot_cmh_adapted_OBO,
             grid_plot_cmh_adapted_OB,
             grid_plot_cmh_adapted_O,
             layout_matrix = layout)

dev.off()

# 4.4 - Scaled and FDR =========================================================
# Pending all plots!!!
# Grid arrange
layout <- matrix(c(1,2,3), ncol = 1, byrow = TRUE)
y_limit_up <- 220

grid_plot_cmh_adapted_OBO <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_obo01_vs_obo20_scaled, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = NULL,
                   x_label = FALSE,
                   y_label = NULL,
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

grid_plot_cmh_adapted_OB <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_ob01_vs_ob20_scaled, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = NULL,
                   x_label = FALSE,
                   y_label = NULL,
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

grid_plot_cmh_adapted_O <-
  GetManhattanPlot(my_dataframe = cmh_pvals,
                   Y = -log10(p.adjust(cmh_pvals$cmh_adapted_o01_vs_o20_scaled, method = "BH")),
                   permutation_pvals = NULL,
                   percentage_significance = TRUE,
                   title = NULL,
                   x_label = FALSE,
                   y_label = NULL,
                   palette = "blue",
                   y_limit_up = y_limit_up,
                   y_limit_down = 0)

# Save grid
png(filename = "results/figures/cmh/adapted_fdr_scaled_O_piled.png",
    width = 1200,
    height = 650)

grid.arrange(grid_plot_cmh_adapted_OBO,
             grid_plot_cmh_adapted_OB,
             grid_plot_cmh_adapted_O,
             layout_matrix = layout)

dev.off()


# 4.4 - QQ Plots ================================================================

n <- nrow(cmh_pvals)

# Create a vector of N values evenly spaces from 1 to 1 / N
N <- sort(-log10(seq(1, n) / n))


p <- sort(-log10(cmh_pvals$cmh_adapted_o01_vs_o20_scaled))

# 5 - PCA Analysis =============================================================
snp_table_shahrestani <- 
  as.data.frame(fread("data/processed/processed_snps_abcd_shahrestani.csv"))

pca_data <- PreparePca(snp_table_shahrestani)

pca_plot_labeled <- PlotPca(pca_data, label = TRUE)
ggsave("results/figures/pca_plot_labeled.png",
       plot = pca_plot_labeled,
       width = 2000,
       height = 2000,
       units = "px")

pca_plot_unlabeled <- PlotPca(pca_data, label = FALSE)
ggsave("results/figures/pca_plot_unlabeled.png",
       plot = pca_plot_unlabeled,
       width = width,
       height = 1600,
       units = "px")

# 6 - Frequency Analysis =======================================================
# Load files
snp_table_shahrestani <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

snp_table_regimes <-
  readRDS("data/processed/processed_snps_abcd_regimes.rds")

cmh_pvals <- readRDS("results/cmh_pvals.rds")

# Check freqs for all SNPs
freq <- GetFreq(snp_table_shahrestani)

# Create a layout with two rows and 5 columns
layout <- matrix(c(1,2,3,4,5,6,7,8,9,10), ncol = 5, byrow = TRUE)

# Plot OBO freqs
obo_rep01_gen01_freq <- 
  ggplot(freq, aes(alt_OBO_rep01_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep02_gen01_freq <-
  ggplot(freq, aes(alt_OBO_rep02_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep03_gen01_freq <-
  ggplot(freq, aes(alt_OBO_rep03_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep04_gen01_freq <-
  ggplot(freq, aes(alt_OBO_rep04_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep05_gen01_freq <-
  ggplot(freq, aes(alt_OBO_rep05_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep01_gen20_freq <-
  ggplot(freq, aes(alt_OBO_rep01_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep02_gen20_freq <-
  ggplot(freq, aes(alt_OBO_rep02_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep03_gen20_freq <-
  ggplot(freq, aes(alt_OBO_rep03_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep04_gen20_freq <-
  ggplot(freq, aes(alt_OBO_rep04_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep05_gen20_freq <-
  ggplot(freq, aes(alt_OBO_rep05_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

png(filename = "results/figures/freq_hist/obo_freq_hist_all.png",
    width = 1800,
    height = 900)

grid.arrange(obo_rep01_gen01_freq,
             obo_rep02_gen01_freq,
             obo_rep03_gen01_freq,
             obo_rep04_gen01_freq,
             obo_rep05_gen01_freq,
             obo_rep01_gen20_freq,
             obo_rep02_gen20_freq,
             obo_rep03_gen20_freq,
             obo_rep04_gen20_freq,
             obo_rep05_gen20_freq,
             layout_matrix = layout)

dev.off()

# Plot OB freqs
OB_rep01_gen01_freq <- 
  ggplot(freq, aes(alt_OB_rep01_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep02_gen01_freq <-
  ggplot(freq, aes(alt_OB_rep02_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep03_gen01_freq <-
  ggplot(freq, aes(alt_OB_rep03_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep04_gen01_freq <-
  ggplot(freq, aes(alt_OB_rep04_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep05_gen01_freq <-
  ggplot(freq, aes(alt_OB_rep05_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep01_gen20_freq <-
  ggplot(freq, aes(alt_OB_rep01_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep02_gen20_freq <-
  ggplot(freq, aes(alt_OB_rep02_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep03_gen20_freq <-
  ggplot(freq, aes(alt_OB_rep03_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep04_gen20_freq <-
  ggplot(freq, aes(alt_OB_rep04_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep05_gen20_freq <-
  ggplot(freq, aes(alt_OB_rep05_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

png(filename = "results/figures/freq_hist/OB_freq_hist_all.png",
    width = 1800,
    height = 900)

grid.arrange(OB_rep01_gen01_freq,
             OB_rep02_gen01_freq,
             OB_rep03_gen01_freq,
             OB_rep04_gen01_freq,
             OB_rep05_gen01_freq,
             OB_rep01_gen20_freq,
             OB_rep02_gen20_freq,
             OB_rep03_gen20_freq,
             OB_rep04_gen20_freq,
             OB_rep05_gen20_freq,
             layout_matrix = layout)

dev.off()

# Now do the same thing for only significant SNPs
# Check freqs for all SNPs
# Add the CMH vals to the snp_table so we can filter them all together
cmh_pvals$ABS_POS <- NULL
cmh_pvals$CHROM <- NULL
snp_table_shahrestani <- cbind(snp_table_shahrestani, cmh_pvals)

threshold <- 1e-100

filtered_snp_table <- snp_table_shahrestani[which(snp_table_shahrestani$cmh_adapted_o01_vs_o20 < threshold),]

freq_significant <- GetFreq(filtered_snp_table)

# Create a layout with two rows and 5 columns
layout <- matrix(c(1,2,3,4,5,6,7,8,9,10), ncol = 5, byrow = TRUE)

# Plot OBO freqs
obo_rep01_gen01_freq_significant <- 
  ggplot(freq_significant, aes(alt_OBO_rep01_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep02_gen01_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep02_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep03_gen01_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep03_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep04_gen01_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep04_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep05_gen01_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep05_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep01_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep01_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep02_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep02_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep03_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep03_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep04_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep04_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

obo_rep05_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OBO_rep05_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

png(filename = "results/figures/freq_hist/obo_freq_significant_hist.png",
    width = 1800,
    height = 900)

grid.arrange(obo_rep01_gen01_freq_significant,
             obo_rep02_gen01_freq_significant,
             obo_rep03_gen01_freq_significant,
             obo_rep04_gen01_freq_significant,
             obo_rep05_gen01_freq_significant,
             obo_rep01_gen20_freq_significant,
             obo_rep02_gen20_freq_significant,
             obo_rep03_gen20_freq_significant,
             obo_rep04_gen20_freq_significant,
             obo_rep05_gen20_freq_significant,
             layout_matrix = layout)

dev.off()

# Plot OB freqs
OB_rep01_gen01_freq_significant <- 
  ggplot(freq_significant, aes(alt_OB_rep01_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep02_gen01_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep02_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep03_gen01_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep03_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep04_gen01_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep04_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep05_gen01_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep05_gen01)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep01_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep01_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep02_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep02_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep03_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep03_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep04_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep04_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

OB_rep05_gen20_freq_significant <-
  ggplot(freq_significant, aes(alt_OB_rep05_gen20)) + 
  geom_histogram(binwidth = 0.01, color = "black") +
  theme_minimal()

png(filename = "results/figures/freq_hist/OB_freq_significant_hist.png",
    width = 1800,
    height = 900)

grid.arrange(OB_rep01_gen01_freq_significant,
             OB_rep02_gen01_freq_significant,
             OB_rep03_gen01_freq_significant,
             OB_rep04_gen01_freq_significant,
             OB_rep05_gen01_freq_significant,
             OB_rep01_gen20_freq_significant,
             OB_rep02_gen20_freq_significant,
             OB_rep03_gen20_freq_significant,
             OB_rep04_gen20_freq_significant,
             OB_rep05_gen20_freq_significant,
             layout_matrix = layout)

dev.off()


# 7 - Allele Trajectory Analysis ===============================================

# Data loading
snp_table_shahrestani <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

cmh_pvals <- readRDS("results/cmh_pvals.rds")


PhaseSnps <- function(freq) {
  
  print("Warning: phasing SNPs put them on a different col order")
  print("Now, we have all gen01 columns and then all gen20 columns")
  
  phased_snp_table_gen01 <- freq[,grep("gen01", colnames(freq))]
  phased_snp_table_gen20 <- freq[,grep("gen20|gen56", colnames(freq))]
  
  marked_for_phasing <- phased_snp_table_gen01 > 0.5
  
  # If marked for phasing, subtract 1 from the value
  phased_snp_table_gen01[marked_for_phasing] <- 
    1 - phased_snp_table_gen01[marked_for_phasing]
  
  phased_snp_table_gen20[marked_for_phasing] <- 
    1 - phased_snp_table_gen20[marked_for_phasing]
  
  phased_snp_table <- cbind(phased_snp_table_gen01, phased_snp_table_gen20)
  
  return(phased_snp_table)
}

cmh_pvals$ABS_POS <- NULL
cmh_pvals$CHROM <- NULL
snp_table_shahrestani <- cbind(snp_table_shahrestani, cmh_pvals)

threshold <- 1e-100

filtered_snp_table <- snp_table_shahrestani[which(snp_table_shahrestani$cmh_adapted_o01_vs_o20 < threshold),]
filtered_snp_table <- dplyr::select(filtered_snp_table, -contains("nB"))

freq_significant <- GetFreq(filtered_snp_table)

phased_freq_significant <- PhaseSnps(freq_significant)

#phased_freq_significant_OBO <- phased_freq_significant[,grep("OBO", colnames(phased_freq_significant))]
#phased_freq_significant_OB <- phased_freq_significant[,grep("OB_", colnames(phased_freq_significant))]

plots <- list()
for (i in 0:9) {
  # Define y axis
  rep_i_01 <- phased_freq_significant[1+i]
  print(colnames(rep_i_01))
  rep_i_20 <- phased_freq_significant[11+i]
  print(colnames(rep_i_20))
  
  colnames(rep_i_01) <- "freq"
  colnames(rep_i_20) <- "freq"
  
  plotting_df <- cbind(rep_i_01, rep_i_20)
  
  # png(filename = paste("results/figures/allele_trajectory/rep", i+1, ".png", sep = ""),
  #     width = 1800,
  #     height = 900)
  
  plot(1,
       type = "n",
       xlab = "", 
       ylab = "",
       xlim = c(1, 20),  
       ylim = c(0, 1),
       xaxt = "n",
       main = paste("rep", i+1),
       cex.axis = 2
  )
  for (row in 1:nrow(plotting_df)) {
    lines(x = c(1,20), 
          y = plotting_df[row,])
  }
  
  p <- recordPlot()
  
  # dev.off()
  
  plots[[i+1]] <- p
  print(i+1)
}


# 7.1 - Delta Statistic ========================================================
# Data loading
snp_table_shahrestani <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

freq <- GetFreq(snp_table_shahrestani)
freq <- dplyr::select(freq, -contains("nB"))

phased_freq <- PhaseSnps(freq)

# Create a delta dataframe, which will calculate the delta in frequency for all the replicates
delta_df <- data.frame(matrix(NA, nrow = nrow(freq), ncol = 10))

gen01 <- phased_freq[,grep("gen01", colnames(phased_freq))]
gen20 <- phased_freq[,grep("gen20|gen56", colnames(phased_freq))]

for (i in 1:10) {
  delta_df[,i] <- gen20[,i] - gen01[,i]
}

colnames(delta_df) <- colnames(freq[1:10])
delta_df <- abs(delta_df)

cmh_pvals <- readRDS("results/cmh_pvals.rds")
cmh_pvals$delta <- rowSums(delta_df)/ncol(delta_df)

y_limit_up <- 220
grid_plot_cmh_adapted_O_scaled_fdr <-
  GetManhattanPlot(
    my_dataframe = cmh_pvals,
    Y = -log10(p.adjust(cmh_pvals$cmh_adapted_o01_vs_o20_scaled, method = "BH")),
    permutation_pvals = NULL,
    percentage_significance = TRUE,
    title = "Adapted CMH test, scaled, FDR corrected: O gen01 vs O gen20, colored for delta > 0.40",
    x_label = TRUE,
    y_label = NULL,
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)


# Lets color the top 10% higher delta values in red
grid_plot_cmh_adapted_O_scaled_fdr <-
  grid_plot_cmh_adapted_O_scaled_fdr +
  aes(color = delta > 0.40) +  # Add color mapping
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  labs(color = "Higher Delta")  # Update legend label


grid_plot_cmh_adapted_O_scaled_fdr

# ==============================================================================

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

# Add the CMH vals to the snp_table so we can filter them all together
cmh_pvals$ABS_POS <- NULL
cmh_pvals$CHROM <- NULL
snp_table_shahrestani_scaled_cmh <- cbind(snp_table_shahrestani_scaled, cmh_pvals)

threshold <- 1e-100

significant_snp_table <- snp_table_shahrestani_scaled_cmh[which(p.adjust(snp_table_shahrestani_scaled_cmh$cmh_adapted_o01_vs_o20_scaled,  method = "fdr") < threshold),]

# freq <- GetFreq(filtered_snp_table)
# 
# # Problem: most significant SNPs start from a fixed frequency all O-type populations
# # If that were true, it would mean that the exact same mutation happened in all experimental populations, which is extremely unlikely
# # This is probably an artifact of our data processing workflow
# 
# # out of the significant snps, make a histogram of the 
# 
# # This will show you that the most significant SNPs start from a fixed frequency
# DiagnoseSnps(snp_table_significant)
# 
# # Let's filter out the SNPs that start from a fixed frequency
# filtered_snp_table <- FilterOutFixedSnps(snp_table_shahrestani)
# 
# # Let's see how a manhattan plot with only these SNPs looks like
# # The perm_pval can't be used here, as the permutation test was run with all SNPs
# # I would need to run a new permutation test with only these 407,678 SNPs to be sure.
# filtered_manhplot <-
#   GetManhattanPlot(
#     my_dataframe = filtered_snp_table,
#     Y = -log10(p.adjust(filtered_snp_table$cmh_adapted_o01_vs_o20_scaled,  method = "fdr")),
#     #permutation_pvals = perm_pvals$o,
#     percentage_significance = TRUE,
#     title = "Filtered Manhattan plot - O gen01 vs O gen20",
#     x_label = TRUE,
#     y_label = "-log10(p-value)",
#     palette = "blue",
#     y_limit_up = 300,
#     y_limit_down = 0
#   )
# 
# # Lets paint the SNPs that are above the threshold in red
# filtered_manhplot <- filtered_manhplot +
#   aes(color = -log10(cmh_adapted_o01_vs_o20) > -log10(threshold)) +  # Add color mapping
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   labs(color = "Above Threshold")  # Update legend label
# 
# filtered_manhplot
# 
# # Let's see the how the most significant SNPs look like now
# filtered_snp_table_significant <-
#   filtered_snp_table[which(filtered_snp_table$cmh_adapted_o01_vs_o20 < threshold),]
# 
# DiagnoseSnps(filtered_snp_table_significant)

# Now we have interesting SNPs to look at. Let's do an enrichment analysis
GO_dataframe <- significant_snp_table[,c("CHROM", "POS", "cmh_adapted_o01_vs_o20_scaled")]
GO_dataframe$cmh_adapted_o01_vs_o20_scaled <- -log10(GO_dataframe$cmh_adapted_o01_vs_o20_scaled)
GO_dataframe$coordinate <- paste(GO_dataframe$CHROM, ":", GO_dataframe$POS, sep = "")



# 8 - Enrichment Analysis ======================================================
# Now we can use biomaRt to get a gene list for those regions.
# We can later use that genelist in a website like GOrilla and see if there is any enrichment
ensembl <- useEnsembl(biomart = "ensembl", dataset = "dmelanogaster_gene_ensembl")


# Get the gene list for each peak
genes <-
  getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = c("chromosome_name", "start", "end"),
    values = list(chromosome_name = GO_dataframe$CHROM, start = GO_dataframe$POS, end = GO_dataframe$POS),
    uniqueRows = TRUE,
    mart = ensembl)

genes_filtered <- genes[!duplicated(genes$external_gene_name),]  
genes_filtered <- genes_filtered[!is.na(genes_filtered$external_gene_name),]
genes_filtered <- genes_filtered[!grepl("df_nrg", genes_filtered$ensembl_gene_id),]

genes_filtered <- genes_filtered[!grepl("^CG", genes_filtered$external_gene_name),]
genes_filtered <- genes_filtered[!grepl("^CR", genes_filtered$external_gene_name),]
genes_filtered <- genes_filtered[!grepl("RNA:", genes_filtered$external_gene_name),]

write.table(genes_filtered$external_gene_name, "results/genes_filtered.txt", quote = FALSE, sep = "\n")

# Using bioconductor tools
BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6") # Genome sequences, not exactly what I want
BiocManager::install("org.Dm.eg.db") # Annotation, that's what I want
BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene") # TxDB (transcription)
BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6") # TxDB (transcription)

#library(Drosophila_melanogaster)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

genes_list <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
genes_list


mycords <- significant_snp_table[c("CHROM", "POS")]
mycords$CHROM <- paste0("chr", mycords$CHROM)
mycords <- 
  mycords %>%
  mutate(chrom=CHROM, start=POS, end=POS) %>%
  dplyr::select(chrom, start, end) %>%
  makeGRangesFromDataFrame()

mycords

subset <- subsetByOverlaps(genes_list, mycords)
genes <- as.data.frame(subset)
genes$symbol <- mapIds(org.Dm.eg.db, keys = genes$gene_id, column = "SYMBOL", keytype = "ENSEMBL")

write.table(unique(genes$symbol), "results/genes_filtered.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n")

as.data.frame(org.Dm.egSYMBOL) %>% head

length(unique(genes$symbol))

# GO Term Analysis
BiocManager::install("clusterProfiler")
library(biomaRt)
library(clusterProfiler)
ensembl <- useEnsembl(biomart = "genes", dataset = "dmelanogaster_gene_ensembl")

go_annotations <-
  biomaRt::getBM(
    attributes = c("ensembl_gene_id", "description", "go_id", "name_1006", "go_linkage_type", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = genes$gene_id,
    mart = ensembl
  )

unwanted_go_evidence <- c("TAS", "NAS", "IC", "ND")

reliable_go_annotations <- go_annotations[!go_annotations$go_linkage_type %in% unwanted_go_evidence,]

go_results <-
  enrichGO(
    gene = go_annotations$ensembl_gene_id,,
    OrgDb = org.Dm.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "fdr",
    qvalueCutoff = 0.05
  )

go_results_df <- as.data.frame(go_results)

print(go_results_df)
