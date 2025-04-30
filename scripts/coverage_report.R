# Coverage Report

snp_table <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

cov <- snp_table[,grep("^N_", colnames(snp_table))]

# Initialize
samples <- colnames(cov)

# Calculate
median_cov <- apply(cov, 2, median, na.rm = TRUE)
mean_cov   <- apply(cov, 2, mean, na.rm = TRUE)
sd_cov     <- apply(cov, 2, sd, na.rm = TRUE)
cv_cov     <- sd_cov / mean_cov

# Combine into a dataframe
coverage_report <- data.frame(
  sample = samples,
  median_cov = median_cov,
  mean_cov = mean_cov,
  sd_cov = sd_cov,
  cv = cv_cov,
  stringsAsFactors = FALSE,
  row.names = NULL
)

# View
print(coverage_report)
