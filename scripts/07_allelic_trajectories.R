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
