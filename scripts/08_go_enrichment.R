# Read snp tables
snp_table_shahrestani <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

snp_table_regimes <-
  readRDS("data/processed/processed_snps_abcd_regimes.rds")

snp_table_shahrestani_scaled <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani_scaled.rds")

snp_table_regimes_scaled <-
  readRDS("data/processed/processed_snps_abcd_regimes_scaled.rds")

cmh_pvals <- readRDS("results/cmh_pvals.rds")

# Add the CMH vals to the snp_table so we can filter them all together
cmh_pvals$ABS_POS <- NULL
cmh_pvals$CHROM <- NULL
snp_table_shahrestani_scaled_cmh <- cbind(snp_table_shahrestani_scaled, cmh_pvals)

threshold <- 1e-100

significant_snp_table <- snp_table_shahrestani_scaled_cmh[which(p.adjust(snp_table_shahrestani_scaled_cmh$cmh_adapted_o01_vs_o20_scaled,  method = "fdr") < threshold),]

# freq <- GetFreq(filtered_snp_table)
# 
# # Problem: most significant SNPs start from a fixed frequency all O-type populations
# # If that were true, it would mean that the exact same mutation happened in all experimental populations, which is extremely unlikely
# # This is probably an artifact of our data processing workflow
# 
# # out of the significant snps, make a histogram of the 
# 
# # This will show you that the most significant SNPs start from a fixed frequency
# DiagnoseSnps(snp_table_significant)
# 
# # Let's filter out the SNPs that start from a fixed frequency
# filtered_snp_table <- FilterOutFixedSnps(snp_table_shahrestani)
# 
# # Let's see how a manhattan plot with only these SNPs looks like
# # The perm_pval can't be used here, as the permutation test was run with all SNPs
# # I would need to run a new permutation test with only these 407,678 SNPs to be sure.
# filtered_manhplot <-
#   GetManhattanPlot(
#     my_dataframe = filtered_snp_table,
#     Y = -log10(p.adjust(filtered_snp_table$cmh_adapted_o01_vs_o20_scaled,  method = "fdr")),
#     #permutation_pvals = perm_pvals$o,
#     percentage_significance = TRUE,
#     title = "Filtered Manhattan plot - O gen01 vs O gen20",
#     x_label = TRUE,
#     y_label = "-log10(p-value)",
#     palette = "blue",
#     y_limit_up = 300,
#     y_limit_down = 0
#   )
# 
# # Lets paint the SNPs that are above the threshold in red
# filtered_manhplot <- filtered_manhplot +
#   aes(color = -log10(cmh_adapted_o01_vs_o20) > -log10(threshold)) +  # Add color mapping
#   scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
#   labs(color = "Above Threshold")  # Update legend label
# 
# filtered_manhplot
# 
# # Let's see the how the most significant SNPs look like now
# filtered_snp_table_significant <-
#   filtered_snp_table[which(filtered_snp_table$cmh_adapted_o01_vs_o20 < threshold),]
# 
# DiagnoseSnps(filtered_snp_table_significant)

# Now we have interesting SNPs to look at. Let's do an enrichment analysis
GO_dataframe <- significant_snp_table[,c("CHROM", "POS", "cmh_adapted_o01_vs_o20_scaled")]
GO_dataframe$cmh_adapted_o01_vs_o20_scaled <- -log10(GO_dataframe$cmh_adapted_o01_vs_o20_scaled)
GO_dataframe$coordinate <- paste(GO_dataframe$CHROM, ":", GO_dataframe$POS, sep = "")



# 8 - Enrichment Analysis ======================================================
# Now we can use biomaRt to get a gene list for those regions.
# We can later use that genelist in a website like GOrilla and see if there is any enrichment
ensembl <- useEnsembl(biomart = "ensembl", dataset = "dmelanogaster_gene_ensembl")


# Get the gene list for each peak
genes <-
  getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = c("chromosome_name", "start", "end"),
    values = list(chromosome_name = GO_dataframe$CHROM, start = GO_dataframe$POS, end = GO_dataframe$POS),
    uniqueRows = TRUE,
    mart = ensembl)

genes_filtered <- genes[!duplicated(genes$external_gene_name),]  
genes_filtered <- genes_filtered[!is.na(genes_filtered$external_gene_name),]
genes_filtered <- genes_filtered[!grepl("df_nrg", genes_filtered$ensembl_gene_id),]

genes_filtered <- genes_filtered[!grepl("^CG", genes_filtered$external_gene_name),]
genes_filtered <- genes_filtered[!grepl("^CR", genes_filtered$external_gene_name),]
genes_filtered <- genes_filtered[!grepl("RNA:", genes_filtered$external_gene_name),]

write.table(genes_filtered$external_gene_name, "results/genes_filtered.txt", quote = FALSE, sep = "\n")

# # Using bioconductor tools
# BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm6") # Genome sequences, not exactly what I want
# BiocManager::install("org.Dm.eg.db") # Annotation, that's what I want
# BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene") # TxDB (transcription)
# BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6") # TxDB (transcription)

#library(Drosophila_melanogaster)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

genes_list <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
genes_list


mycords <- significant_snp_table[c("CHROM", "POS")]
mycords$CHROM <- paste0("chr", mycords$CHROM)
mycords <- 
  mycords %>%
  mutate(chrom=CHROM, start=POS, end=POS) %>%
  dplyr::select(chrom, start, end) %>%
  makeGRangesFromDataFrame()

mycords

subset <- subsetByOverlaps(genes_list, mycords)
genes <- as.data.frame(subset)
genes$symbol <- mapIds(org.Dm.eg.db, keys = genes$gene_id, column = "SYMBOL", keytype = "ENSEMBL")

write.table(unique(genes$gene_id), "results/genes_filtered_id.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\n")

as.data.frame(org.Dm.egSYMBOL) %>% head

length(unique(genes$symbol))

# GO Term Analysis
# BiocManager::install("clusterProfiler")
library(biomaRt)
library(clusterProfiler)
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
    gene = go_annotations$ensembl_gene_id,,
    OrgDb = org.Dm.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "fdr",
    qvalueCutoff = 0.05
  )

go_results_df <- as.data.frame(go_results)

simple_go_results <- simplify(go_results)
simple_go_results_df <- as.data.frame(simple_go_results)


print(go_results_df)
dotplot(go_results, showCategory = 20)
dotplot(simple_go_results, showCategory = 20)
