# B annotation
# Clean environment
source("scripts/utils.R")
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(clusterProfiler)
library(biomaRt)

snp_table <- readRDS("data/processed/processed_snps_shahrestani_scaled.rds") # Note that this is shahrestani now
all_cmh <- readRDS("results/all_cmh.rds")

all_cmh <- all_cmh[,6:ncol(all_cmh)]
snp_table_cmh <- cbind(snp_table, all_cmh)

# Save RAM
rm(all_cmh)
rm(snp_table)

column <- snp_table_cmh$adapted_nbo56_vs_nbo01_fdr

threshold <- 1e-100
#threshold <- quantile(column, 0.001)

significant_snp_table <- snp_table_cmh[which(column < threshold),]

mycords <- significant_snp_table[c("CHROM", "POS")]
mycords$CHROM <- paste0("chr", mycords$CHROM)

mycords <- 
  mycords %>%
  mutate(chrom=CHROM, start=POS, end=POS) %>%
  dplyr::select(chrom, start, end) %>%
  makeGRangesFromDataFrame()

subset <- subsetByOverlaps(genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene), mycords)
genes <- as.data.frame(subset) # 0 genes, nothing to do here