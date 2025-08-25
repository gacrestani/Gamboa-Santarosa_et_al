long_genes <- long_genes$V1

intersect(long_genes, genes$gene_id)




library(readxl)
carnes <- read_excel("~/Downloads/pone.0138569.s014.xlsx")

intersect(genes$gene_id, carnes$`FlyBase ID`)


fabian <- read_excel("~/Downloads/evl389-sup-0002-tables1.xls", sheet = "candidate SNPs", skip = 1)

intersect_fabian <- intersect(genes$gene_id, unique(fabian$`Flybase ID`))

genes_fabian <- genes[genes$gene_id %in% intersect_fabian, ]
genes_fabian$symbol

cat(sort(genes_fabian$symbol), sep = ", ")

intersect(intersect_fabian, carnes$`FlyBase ID`)

intersect(intersect(fabian$`Flybase ID`, carnes$`FlyBase ID`), genes$gene_id)


carnes[carnes$`FlyBase ID` %in% intersect(fabian$`Flybase ID`, carnes$`FlyBase ID`), ]
