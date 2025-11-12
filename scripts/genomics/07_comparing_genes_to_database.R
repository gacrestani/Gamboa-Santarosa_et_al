source("scripts/utils.R")

library(org.Dm.eg.db)
library(GO.db)

get_all_genes_from_GO <- function(parent_go_id) {
  all_bp_offspring <- as.list(GOBPOFFSPRING)
  child_go_id <- all_bp_offspring[[parent_go_id]]
  
  all_go_ids <- c(parent_go_id, child_go_id)
  
  # Query dmel database
  gene_annotations <- AnnotationDbi::select(org.Dm.eg.db,
                                            keys = all_go_ids,
                                            columns = c("FLYBASE", "SYMBOL", "GENENAME"),
                                            keytype = "GO")
  
  return(gene_annotations)
}

genes <- read.csv("results/significant_genes.csv")

# IMMUNE RESPONSE
lifespan_response_genes <- get_all_genes_from_GO(parent_go_id = "GO:0008340") # determination of adult lifespan
immune_response_genes <- get_all_genes_from_GO(parent_go_id = "GO:0006955") # immune response


lifespan_related_genes <- genes[genes$gene_id %in% lifespan_response_genes$FLYBASE, ]
immune_related_genes <- genes[genes$gene_id %in% immune_response_genes$FLYBASE, ]