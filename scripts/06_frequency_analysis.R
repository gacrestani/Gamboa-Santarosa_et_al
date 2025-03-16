# frequency_analysis.R
# Analyses frequencies of my samples, looking at deltas and other metrics
#
# inputs: snp_table.rds (x4), cmh_pvals.rds, perm_pvals.csv
# outputs: manhattan plots, in accordance with the following format:
#
# THIS NEEDS CLEANING!

source("scripts/functions.R")

# 0 Load files  ----------------------------------------------------------------
snp_table_shahrestani <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

snp_table_regimes <-
  readRDS("data/processed/processed_snps_abcd_regimes.rds")

cmh_pvals <- readRDS("results/cmh_pvals.rds")

# Parameters
width      <- 720
height     <- 480 

treatments <- c("OBO", "OB", "nBO", "nB")
gen2 <- c("20", "20", "56", "56")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
