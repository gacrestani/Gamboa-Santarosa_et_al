# cmh_tests.R
# Runs the CMH tests on the SNP tables.
# It runs the classic and adapted CMH test on scaled and untransformed snp tables.
#
# inputs: snp_table.rds (x4)
# outputs: cmh_pvals.rds - a dataframe with the p-values for each test.

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

# Initiates an empty dataframe with the same number of rows as our snp tables
cmh_pvals <-
  data.frame(matrix(NA,
                    nrow = nrow(snp_table_shahrestani),
                    ncol = 0))

cmh_pvals$ABS_POS <- snp_table_shahrestani$ABS_POS
cmh_pvals$CHROM <- snp_table_shahrestani$CHROM

treatments <- c("OBO", "OB", "nBO", "nB", "O", "B")
gen2 <- c("20", "20", "56", "56", "20", "56")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 1 Classic CMH ----------------------------------------------------------------
# NOTE: These functions take a while to run (even running in parallel)
for (i in 1:length(treatments)) {
  
  # Define if using regimes or shahrestani snp table
  ifelse(
    treatments[i] == "O" | treatments[i] == "B",
    snp_table_iter <- snp_table_regimes,
    snp_table_iter <- snp_table_shahrestani
  )
  
  # Run test
  cmh_pvals[[paste0("cmh_classic_", treatments[i])]] <-
    ClassicalCmhTest(
      snp_table = snp_table_iter,
      treatment1 = treatments[i])
  
  # Print test name completed
  print(paste("Classic CMH test done:", treatments[i]))
  
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 2 Adapted CMH ----------------------------------------------------------------
# These are faster than the classical tests
for (i in 1:length(treatments)) {
  
  # Define if using regimes or shahrestani snp table
  ifelse(
    treatments[i] == "O" | treatments[i] == "B",
    snp_table_iter <- snp_table_regimes,
    snp_table_iter <- snp_table_shahrestani
  )
  
  # Run test
  cmh_pvals[[paste0("cmh_adapted_", treatments[i])]] <-
    AdaptedCmhTest(
      snp_table = snp_table_iter,
      treatment1 = treatments[i])
  
  # Print test name completed
  print(paste("Classic CMH test done:", treatments[i]))
  
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 3 Classic CMH ----------------------------------------------------------------
# Takes a while to run
for (i in 1:length(treatments)) {
  
  # Define if using regimes or shahrestani snp table
  ifelse(
    treatments[i] == "O" | treatments[i] == "B",
    snp_table_iter <- snp_table_regimes_scaled,
    snp_table_iter <- snp_table_shahrestani_scaled
  )
  
  # Run test
  cmh_pvals[[paste0("cmh_classic_scaled_", treatments[i])]] <-
    ClassicalCmhTest(
      snp_table = snp_table_iter,
      treatment1 = treatments[i])
  
  # Print test name completed
  print(paste("Classic CMH test done:", treatments[i]))

}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 4 Adapted CMH ----------------------------------------------------------------
for (i in 1:length(treatments)) {
  
  # Define if using regimes or shahrestani snp table
  ifelse(
    treatments[i] == "O" | treatments[i] == "B",
    snp_table_iter <- snp_table_regimes_scaled,
    snp_table_iter <- snp_table_shahrestani_scaled
  )
  
  # Run test
  cmh_pvals[[paste0("cmh_adapted_scaled_", treatments[i])]] <-
    AdaptedCmhTest(
      snp_table = snp_table_iter,
      treatment1 = treatments[i])
  
  # Print test name completed
  print(paste("Classic CMH test done:", treatments[i]))
  
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 5 FDR correct everything -----------------------------------------------------

fdr_cmh <- as.data.frame(lapply(cmh_pvals[-c(1,2)], function(p) p.adjust(p, method = "fdr")))
names(fdr_cmh) <- sub("(cmh_[^_]+_)", "\\1fdr_", names(cmh_pvals)[-c(1,2)])

all_cmh <- cbind(cmh_pvals, fdr_cmh)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 5 Save results ---------------------------------------------------------------
saveRDS(all_cmh, "results/cmh_pvals.rds")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~