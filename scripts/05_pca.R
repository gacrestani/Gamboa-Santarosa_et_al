source("scripts/functions.R")

# Hannah's stuff
snp_table <- read.table("/nfs3/IB/Burke_Lab/2025_FLU/2025_02_21_pipeline_run01/results/filtered_snps.txt", header = TRUE)
snp_table[snp_table == "."] <- 0
snp_table[, 6:ncol(snp_table)] <- lapply(snp_table[, 6:ncol(snp_table)], as.numeric)
snp_table <- FilterMinAndMaxCov(snp_table, min_cov = 5, max_cov = 500)

pca <- PreparePca(snp_table)
PlotPca(pca, label = FALSE)

