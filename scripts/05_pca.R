source("scripts/functions.R")

pca <- PreparePca(snp_table_shahrestani)
PlotPca(pca, label = FALSE)