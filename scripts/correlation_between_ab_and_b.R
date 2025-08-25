source("scripts/functions.R")

snp_table_a <- ReadSnpTable(path = "data/snp_tables/filtered_snps_fly02a.txt")
snp_table_b <- ReadSnpTable(path = "data/snp_tables/filtered_snps_fly02b.txt")
snp_table_ab <- ReadSnpTable(path = "data/snp_tables/filtered_snps_fly02ab.txt")

freq_a <- GetFreq(snp_table_a)
freq_b <- GetFreq(snp_table_b)
freq_ab <- GetFreq(snp_table_ab)

freq_a$id <- paste0(snp_table_a$CHROM, ":", snp_table_a$POS)
freq_b$id <- paste0(snp_table_b$CHROM, ":", snp_table_b$POS)
freq_ab$id <- paste0(snp_table_ab$CHROM, ":", snp_table_ab$POS)

freq_a <- freq_a[complete.cases(freq_a),]
freq_b <- freq_b[complete.cases(freq_b),]
freq_ab <- freq_ab[complete.cases(freq_ab),]

#common_snps <- intersect(freq_ab$id, freq_b$id)
common_snps <- intersect(freq_a$id, freq_b$id)

freq_a <- freq_a[freq_a$id %in% common_snps,]
freq_b <- freq_b[freq_b$id %in% common_snps,]
freq_ab <- freq_ab[freq_ab$id %in% common_snps,]

freq_a$id <- NULL
freq_b$id <- NULL
freq_ab$id <- NULL

#cor_list <- mapply(function(x, y) cor(x, y), freq_ab, freq_b)
cor_list <- mapply(function(x, y) cor(x, y), freq_a, freq_b)

