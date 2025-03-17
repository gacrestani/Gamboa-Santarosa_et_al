source("scripts/functions.R")

snp_table_shahrestani <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

raw_snp_table <- as.data.frame(fread("data/snp_tables/filtered_snps_abcd.txt"))
raw_snp_table <- AddAbsPosToSnpTable(raw_snp_table)
cov_raw <- GetCov(raw_snp_table)

cov <- GetCov(snp_table_shahrestani)

plot(cov$N_nBO_rep01_gen01 ~ snp_table_shahrestani$ABS_POS, type = "l")
plot(cov_raw$N_CBO_rep01_gen01 ~ raw_snp_table$ABS_POS, type = "l")
