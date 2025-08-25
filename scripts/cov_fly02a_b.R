# Coverage Report



snp_table_a <-
  ReadSnpTable(path = "/var/home/crestang/Downloads/filtered_snps_FLY02a.txt")

snp_table_b <-
  ReadSnpTable(path = "/var/home/crestang/Downloads/filtered_snps_FLY02b.txt")


cov_a <- snp_table_a[,grep("^N_", colnames(snp_table_a))]
cov_b <- snp_table_b[,grep("^N_", colnames(snp_table_b))]

cov <- cov_a
cov <- cov_b

# Initialize
samples <- colnames(cov)

# Calculate
median_cov <- apply(cov, 2, median, na.rm = TRUE)
mean_cov   <- apply(cov, 2, mean, na.rm = TRUE)
sd_cov     <- apply(cov, 2, sd, na.rm = TRUE)
cv_cov     <- sd_cov / mean_cov

# Combine into a dataframe
coverage_report2 <- data.frame(
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


# Normalized coverage
norm_cov <- sweep(cov, 2, colMeans(cov, na.rm = TRUE), FUN = "/")


# Plot each column
for (i in colnames(norm_cov)) {
  png(paste0("results/figures/coverage/normalized_coverage_", colnames(norm_cov[i]), ".png"))
  plot(norm_cov[[i]], main = colnames(norm_cov[i]))
  abline(v = which(snp_table$POS == 3704168), col = "red")
  abline(v = which(snp_table$POS == 4852295), col = "red")
  abline(v = which(snp_table$POS == 25211163), col = "red")
  dev.off()
}


which(snp_table$POS == 3702504) - which(snp_table$POS == 3704733)
which(snp_table$POS == 4853176) - which(snp_table$POS == 4849264)
which(snp_table$POS == 25219680) - which(snp_table$POS == 25209917)
