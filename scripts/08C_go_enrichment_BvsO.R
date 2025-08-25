source("scripts/functions.R")

snp_table_scaled <- 
  readRDS("data/processed/processed_snps_abcd_regimes_scaled.rds")

cmh_pvals <- readRDS("results/cmh_pvals.rds")

cmh_pvals$cmh_classic_o20_vs_b55_scaled <- ClassicalCmhTest(
  snp_table_scaled,
  treatment1 = "O",
  gen1 = "20",
  treatment2 = "B",
  gen2 = "56"
)

cmh_pvals$cmh_classic_fdr_o20_vs_b55_scaled <- p.adjust(cmh_pvals$cmh_classic_o20_vs_b55_scaled, method = "fdr")

cmh_pvals$ABS_POS <- NULL
cmh_pvals$CHROM <- NULL
snp_table_scaled_cmh <- cbind(snp_table_scaled, cmh_pvals)

threshold <- 1e-100

pvals <- cmh_pvals$cmh_classic_fdr_o20_vs_b55_scaled

significant_snp_table <- snp_table_scaled_cmh[which(pvals < threshold),]

library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

genes_list <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
genes_list

mycords <- significant_snp_table[c("CHROM", "POS")]
mycords$CHROM <- paste0("chr", mycords$CHROM)

#colnames(mycords) <- c("chr", "pos")

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