source("scripts/utils.R")
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(VariantAnnotation)
library(clusterProfiler)
library(biomaRt)

snp_table <- readRDS("data/processed/processed_snps_regimes_scaled.rds")
all_cmh <- readRDS("results/all_cmh.rds")

# Add the CMH vals to the snp_table so we can filter them all together

all_cmh <- all_cmh[,6:ncol(all_cmh)]
snp_table_cmh <- cbind(snp_table, all_cmh)

# Save RAM
rm(all_cmh)
rm(snp_table)

# Choose which comparison you will focus on
column <- snp_table_cmh$adapted_o20_vs_o01_fdr

threshold <- 1e-100
#threshold <- quantile(column, 0.001)

significant_snp_table <- snp_table_cmh[which(column < threshold),]

mycords_df <- significant_snp_table[c("CHROM", "POS", "REF", "ALT")]
mycords_df$CHROM <- paste0("chr", mycords_df$CHROM)

mycords_df <- 
  mycords_df %>%
  mutate(chrom=CHROM, start=POS, end=POS, ref=REF, alt=ALT) %>%
  dplyr::select(chrom, start, end, ref, alt)

rownames(mycords_df) <- paste0(mycords_df$chrom, ":", mycords_df$start, "-", mycords_df$end)

sig_snps <- GRanges(seqnames = mycords_df$chrom,
                    ranges = IRanges(start = mycords_df$start, end = mycords_df$end),
                    strand = "*",
                    ref = mycords_df$ref,
                    alt = mycords_df$alt)


subset <- subsetByOverlaps(genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene), sig_snps)
genes <- as.data.frame(subset)
genes$symbol <- mapIds(org.Dm.eg.db, keys = genes$gene_id, column = "SYMBOL", keytype = "ENSEMBL")

# Annotate significant SNPs
loc_annotation <- locateVariants(query = sig_snps,
                                 subject = TxDb.Dmelanogaster.UCSC.dm6.ensGene,
                                 region = AllVariants())
annotation <- as.data.frame(loc_annotation)

list_columns <- sapply(annotation, is.list)
annotation[list_columns] <- lapply(annotation[list_columns], function(col) {
  sapply(col, paste, collapse = ";")
})




# Done capturing genes. Now let's do the GO enrichment analysis
ensembl <- useEnsembl(biomart = "genes", dataset = "dmelanogaster_gene_ensembl")

# Get GO terms associated with my genes
go_annotations <-
  biomaRt::getBM(
    attributes = c("ensembl_gene_id", "description", "go_id", "name_1006", "go_linkage_type", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = genes$gene_id,
    mart = ensembl
  )

unwanted_go_evidence <- c("TAS", "NAS", "IC", "ND")
reliable_go_annotations <- go_annotations[!go_annotations$go_linkage_type %in% unwanted_go_evidence,]

aggregated_go <- reliable_go_annotations %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    go_ids = paste(unique(go_id), collapse = "; "),
    go_names = paste(unique(name_1006), collapse = "; "),
    descriptions = paste(unique(description), collapse = "; "),
    external_gene_names = paste(unique(external_gene_name), collapse = "; ")
  )

gene_list_result <- left_join(genes, aggregated_go, by = c("gene_id" = "ensembl_gene_id"))

# Perform GO enrichment analysis
go_results <-
  enrichGO(
    gene = unique(reliable_go_annotations$ensembl_gene_id),
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
dotplot <- dotplot(simple_go_results, showCategory = 21) + ggtitle("GO Enrichment Analysis - non-redundant terms")

# Save outputs

# Gene list
write.table(gene_list_result, file = "results/significant_genes.csv", sep = "\t", row.names = FALSE, quote = FALSE)

# GO results
write.table(go_results_df, file = "results/go_enrichment_results.csv", sep = "\t", row.names = FALSE, quote = FALSE)

# Annotation
write.table(annotation, file = "results/significant_snp_annotations.csv", sep = "\t", row.names = FALSE, quote = FALSE)

# Save dotplot
ggsave(
  "results/figures/figure5.jpeg",
  dotplot,
  width = 8,
  height = 11,
  units = "in",
  dpi = 300
)
