source("scripts/utils.R")

genes <- read.csv("results/significant_genes.csv")

carnes <- readxl::read_excel("data/cited_papers/pone.0138569.s014.xlsx")
intersection_with_carnes <- genes[genes$gene_id %in% carnes$`FlyBase ID`, ]

fabian <- readxl::read_excel("data/cited_papers/evl389-sup-0002-tables1.xls", sheet = "candidate SNPs", skip = 1)
fabian <- fabian[-1,]

intersection_with_fabian <- genes[genes$gene_id %in% fabian$`Flybase ID`, ]

# Three way intersection
fabian_and_carnes <- carnes[carnes$`FlyBase ID` %in% fabian$`Flybase ID`, ]

three_intersection <- genes[genes$gene_id %in% fabian_and_carnes$`FlyBase ID`, ]