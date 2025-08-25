source("scripts/functions.R")

snp_table <- # Use this for B1-10 gen 01 vs B1-10 gen 55
  readRDS("data/processed/processed_snps_abcd_regimes_scaled.rds")

snp_table <- # Use this for nBO1-5 gen01 vs nBO1-5 gen 55
  readRDS("data/processed/processed_snps_abcd_shahrestani_scaled.rds")

cmh_pvals <- readRDS("results/cmh_pvals.rds")
cmh_pvals$ABS_POS <- NULL
cmh_pvals$CHROM <- NULL
snp_table_cmh <- cbind(snp_table, cmh_pvals)

peak_identification <- snp_table_cmh %>% dplyr::select(CHROM, POS, ABS_POS, cmh_adapted_b01_vs_b56_scaled)
peak_identification$log10 <- -log10(peak_identification$cmh_adapted_b01_vs_b56_scaled)

# B1-10 peaks:
# Peak1: start = 2R:3,702,504 max: 2R:3,704,168 (64) end: 2R:3,704,733
# Peak2: start = 2R:4,849,264 max: 2R:4,852,295 (91) end: 2R:4,853,176
# Peak3: start = 3L:25,209,917 max: 3L:25,211,163 (91) end: 3L:25,219,680

column <- snp_table_cmh$cmh_adapted_nbo01_vs_nbo56_scaled
column <- snp_table_cmh$cmh_adapted_b01_vs_b56_scaled

significant_snp_table2 <- snp_table_cmh[which(p.adjust(column,  method = "fdr") < threshold),]

library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

genes_list <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

mycords <- significant_snp_table[c("CHROM", "POS")]
mycords$CHROM <- paste0("chr", mycords$CHROM)

mycords <- 
  mycords %>%
  mutate(chrom=CHROM, start=POS, end=POS) %>%
  dplyr::select(chrom, start, end) %>%
  makeGRangesFromDataFrame()

subset <- subsetByOverlaps(genes_list, mycords)
genes <- as.data.frame(subset)
genes$symbol <- mapIds(org.Dm.eg.db, keys = genes$gene_id, column = "SYMBOL", keytype = "ENSEMBL")


library(clusterProfiler)
library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "dmelanogaster_gene_ensembl")

go_annotations <-
  biomaRt::getBM(
    attributes = c("ensembl_gene_id", "description", "go_id", "name_1006", "go_linkage_type", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = genes$gene_id,
    mart = ensembl
  )

unwanted_go_evidence <- c("TAS", "NAS", "IC", "ND")

reliable_go_annotations <- go_annotations[!go_annotations$go_linkage_type %in% unwanted_go_evidence,]

go_results <-
  enrichGO(
    gene = reliable_go_annotations$ensembl_gene_id,
    OrgDb = org.Dm.eg.db,
    keyType = "ENSEMBL",
    ont = "ALL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "fdr",
    qvalueCutoff = 0.05
  )

go_results_df <- as.data.frame(go_results)

simple_go_results <- simplify(go_results)
simple_go_results_df <- as.data.frame(simple_go_results)










