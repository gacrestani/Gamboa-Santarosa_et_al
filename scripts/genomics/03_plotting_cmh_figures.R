# plotting_cmh_figures.R

source("scripts/utils.R")

GetManhattanPlot <-
  function(
    my_dataframe,
    Y,
    percentage_significance = NULL,
    title = "Manhattan Plot",
    x_label = TRUE,
    y_label = "-log10(p)",
    palette = "blue",
    y_limit_up = 100,
    y_limit_down = 0) {
    
    # This huge function plots a good manhattan plot considering all the settings I have on the arguments.
    # Some of them have been tweaked over the time to reach what I think is an optimal graph for my use case
    axis_set <-
      c("2L" = (max(my_dataframe[my_dataframe$CHROM == "2L", ]$ABS_POS) + 
                  min(my_dataframe[my_dataframe$CHROM == "2L", ]$ABS_POS))/2,
        "2R" = (max(my_dataframe[my_dataframe$CHROM == "2R", ]$ABS_POS) + 
                  min(my_dataframe[my_dataframe$CHROM == "2R", ]$ABS_POS))/2,
        "3L" = (max(my_dataframe[my_dataframe$CHROM == "3L", ]$ABS_POS) + 
                  min(my_dataframe[my_dataframe$CHROM == "3L", ]$ABS_POS))/2,
        "3R" = (max(my_dataframe[my_dataframe$CHROM == "3R", ]$ABS_POS) + 
                  min(my_dataframe[my_dataframe$CHROM == "3R", ]$ABS_POS))/2,
        "X"  = (max(my_dataframe[my_dataframe$CHROM == "X",  ]$ABS_POS) + 
                  min(my_dataframe[my_dataframe$CHROM == "X",  ]$ABS_POS))/2)
    
    p <- ggplot(my_dataframe,
                aes(x = ABS_POS, 
                    y = Y,
                    color = as.factor(CHROM))) +
      geom_point(alpha = 0.75) +
      scale_x_continuous(labels = names(axis_set), breaks = axis_set) +
      {if (!x_label) scale_x_continuous(labels = NULL, breaks = axis_set)} +
      ylim(y_limit_down, y_limit_up) +
      {if (palette == "blue")
        scale_color_manual(values = rep(c("steelblue1", "steelblue4"), unique(length(axis_set))))} +
      {if (palette == "red")
        scale_color_manual(values = rep(c("firebrick1", "firebrick4"), unique(length(axis_set))))} +
      {if (palette == "grey")
        scale_color_manual(values = rep(c("black", "grey"), unique(length(axis_set))))} +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = NULL, y = y_label) +
      theme_bw() +
      theme(legend.position="none") +
      ggtitle(title)
    
    if (percentage_significance) {
      threshold <- 100 # or -log10(1e-100)
      p <- p +
        geom_hline(yintercept = threshold, col = "red")
    }
    
    return(p)
  }

all_cmh <- readRDS("results/all_cmh.rds")

# Parameters
y_limit_up <- 220 # manhattan plot y-axis upper limit




# Adapted O
adapted_O <-
  GetManhattanPlot(
    my_dataframe = all_cmh,
    Y = -log10(all_cmh$adapted_o20_vs_o01_fdr),
    percentage_significance = TRUE,
    title = expression(paste("O"["1-10"], " generation 20 vs. O"["1-10"], " generation 1")),
    x_label = FALSE,
    y_label = " ",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

adapted_OBO <-
  GetManhattanPlot(
    my_dataframe = all_cmh,
    Y = -log10(all_cmh$adapted_obo20_vs_obo01_fdr),
    percentage_significance = TRUE,
    title = expression(paste("OBO"["1-5"], " generation 20 vs. OBO"["1-5"], " generation 1")),
    x_label = FALSE,
    y_label = expression(paste("-log"["10"], "(FDR-adjusted p-value)")),
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

adapted_OB <-
  GetManhattanPlot(
    my_dataframe = all_cmh,
    Y = -log10(all_cmh$adapted_ob20_vs_ob01_fdr),
    percentage_significance = TRUE,
    title = expression(paste("OB"["1-5"], " generation 20 vs. OB"["1-5"], " generation 1")),
    x_label = TRUE,
    y_label = " ",
    palette = "blue",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

adapted_o_grid <- ggarrange(
  adapted_O,
  adapted_OBO,
  adapted_OB,
  labels = c("A", "B", "C"),
  ncol = 1,
  nrow = 3
)

ggsave(
  "results/adapted_o.png",
  adapted_o_grid,
  width = 11,
  height = 9,
  units = "in",
  dpi = 600
)



# Adapted B
adapted_B <-
  GetManhattanPlot(
    my_dataframe = all_cmh,
    Y = -log10(all_cmh$adapted_b56_vs_b01_fdr),
    percentage_significance = TRUE,
    title = expression(paste("B"["1-10"], " generation 56 vs. B"["1-10"], " generation 1")),
    x_label = FALSE,
    y_label = " ",
    palette = "red",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

adapted_nBO <-
  GetManhattanPlot(
    my_dataframe = all_cmh,
    Y = -log10(all_cmh$adapted_nbo56_vs_nbo01_fdr),
    percentage_significance = TRUE,
    title = expression(paste("nBO"["1-5"], " generation 56 vs. nBO"["1-5"], " generation 1")),
    x_label = FALSE,
    y_label = expression(paste("-log"["10"], "(FDR-adjusted p-value)")),
    palette = "red",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

adapted_nB <-
  GetManhattanPlot(
    my_dataframe = all_cmh,
    Y = -log10(all_cmh$adapted_nb56_vs_nb01_fdr),
    percentage_significance = TRUE,
    title = expression(paste("nB"["1-5"], " generation 56 vs. nB"["1-5"], " generation 1")),
    x_label = TRUE,
    y_label = " ",
    palette = "red",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

adapted_b_grid <- ggarrange(
  adapted_B,
  adapted_nBO,
  adapted_nB,
  labels = c("A", "B", "C"),
  ncol = 1,
  nrow = 3
)

ggsave(
  "results/adapted_b.png",
  adapted_b_grid,
  width = 11,
  height = 9,
  units = "in",
  dpi = 600
)



# Classical B vs O
classical_OvsB <-
  GetManhattanPlot(
    my_dataframe = all_cmh,
    Y = -log10(all_cmh$classical_o20_vs_b56_fdr),
    percentage_significance = TRUE,
    title = expression(paste("Classical CMH - O"["1-10"], " generation 20 vs. B"["1-10"], " generation 56")),
    x_label = FALSE,
    y_label = " ",
    palette = "grey",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

classical_OBOvsnBO <-
  GetManhattanPlot(
    my_dataframe = all_cmh,
    Y = -log10(all_cmh$classical_obo20_vs_nbo56_fdr),
    percentage_significance = TRUE,
    title = expression(paste("Classical CMH - OBO"["1-5"], " generation 20 vs. nBO"["1-5"], " generation 56")),
    x_label = FALSE,
    y_label = expression(paste("-log"["10"], "(FDR-adjusted p-value)")),
    palette = "grey",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

classical_OBvsnB <-
  GetManhattanPlot(
    my_dataframe = all_cmh,
    Y = -log10(all_cmh$classical_ob20_vs_nb56_fdr),
    percentage_significance = TRUE,
    title = expression(paste("Classical CMH - OB"["1-5"], " generation 20 vs. nB"["1-5"], " generation 56")),
    x_label = TRUE,
    y_label = " ",
    palette = "grey",
    y_limit_up = y_limit_up,
    y_limit_down = 0)

classical_o_vs_b <- ggarrange(
  classical_OvsB,
  classical_OBOvsnBO,
  classical_OBvsnB,
  labels = c("A", "B", "C"),
  ncol = 1,
  nrow = 3
)

ggsave(
  "results/classical_o_vs_b.png",
  classical_o_vs_b,
  width = 11,
  height = 9,
  units = "in",
  dpi = 600
)