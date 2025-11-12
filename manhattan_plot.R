GetManhattanPlot <-
  function(
    my_dataframe,
    Y,
    permutation_pvals = NULL,
    percentage_significance = NULL,
    title = "Manhattan Plot",
    x_label = TRUE,
    y_label = "-log10(p)",
    palette = "blue",
    y_limit_up = 220,
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
      ylim(y_limit_down, 220) +
      scale_color_manual(values = rep(c("steelblue1", "steelblue4"), unique(length(axis_set)))) +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = NULL, y = y_label) +
      theme_minimal() +
      theme(legend.position="none") +
      ggtitle(title)
    
    if (percentage_significance) {
      threshold <- 100 # or -log10(1e-100)
      p <- p +
        geom_hline(yintercept = threshold, col = "red")
    }
    
    return(p)
  }


# ------------------------------------------------------------------------------

o <- GetManhattanPlot(
  my_dataframe = cmh_pvals,
  Y = -log10(cmh_pvals$cmh_adapted_o01_vs_o20),
  #permutation_pvals = perm_pvals$o,
  percentage_significance = TRUE,
  title = "O",
  x_label = TRUE,
  y_label = "",
  palette = "blue",
  y_limit_up = 220,
  y_limit_down = 0)

o


sample_size <- floor(0.001 * nrow(cmh_pvals))
sample_indices <- sample(nrow(cmh_pvals), size = sample_size, replace = FALSE)
red_cmh_pvals <- cmh_pvals[sample_indices,]
write.csv(red_cmh_pvals, "~/Analyses/cmh_pvals.csv")
