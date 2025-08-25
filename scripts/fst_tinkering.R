source("scripts/functions.R")

# O 1-10 gen 20 compared to O 1-10 gen 01
snp_table <- # Regimes scaled
  readRDS("data/processed/processed_snps_abcd_regimes.rds")

freqs <- GetFreq(snp_table)
o_freqs <- freqs[,grep("^alt_O_", colnames(freqs))]

p1 <- o_freqs[,grep("*_gen01", colnames(o_freqs))]
p2 <- o_freqs[,grep("*_gen20", colnames(o_freqs))]

q1 <- 1 - p1
q2 <- 1 - p2

pt <- (p1 + p2) / 2
qt <- (q1 + q2) / 2

Ht <- 2 * pt * qt

Hs1 <- 2 * p1 * q1
Hs2 <- 2 * p2 * q2

Hs <- (Hs1 + Hs2) / 2
Fst <- (Ht - Hs) / Ht


Fst <- rowSums(Fst)/ncol(Fst)
Ht <- rowSums(Ht)/ncol(Ht)

snp_table$Fst <- Fst
snp_table$Ht <- Ht


# Fst - O20
o_freqs_20 <- o_freqs[,grep("*_gen20", colnames(o_freqs))]

pn <- o_freqs_20
qn <- 1 - pn

pt <- rowSums(pn) / ncol(pn)
qt <- rowSums(qn) / ncol(qn)

Ht <- 2 * pt * qt

Hsn <-  2 * pn * qn
Hs <- rowSums(Hsn) / ncol(Hsn)

Fst <- (Ht - Hs) / Ht


# Fst - OBO1-5 gen20
snp_table <-
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

freqs <- GetFreq(snp_table)

o_freqs <- freqs[,grep("^alt_OBO_", colnames(freqs))]
o_freqs_20 <- o_freqs[,grep("*_gen20", colnames(o_freqs))]

pn <- o_freqs_20
qn <- 1 - pn

pt <- rowSums(pn) / ncol(pn)
qt <- rowSums(qn) / ncol(qn)

Ht <- 2 * pt * qt

Hsn <-  2 * pn * qn
Hs <- rowSums(Hsn) / ncol(Hsn)

Fst <- (Ht - Hs) / Ht

Fst[is.infinite(Fst)] <- 0
Fst[is.nan(Fst)] <- 0

summary(Fst)
hist(Fst)



# Fst - OB1-5 gen20
snp_table <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")

freqs <- GetFreq(snp_table)

o_freqs <- freqs[,grep("^alt_OB_", colnames(freqs))]
o_freqs_20 <- o_freqs[,grep("*_gen20", colnames(o_freqs))]

GetFst(o_freq_20)





GetFst <- function(freqs) {
  
  pn <- o_freqs_20
  qn <- 1 - pn
  
  pt <- rowSums(pn) / ncol(pn)
  qt <- rowSums(qn) / ncol(qn)
  
  Ht <- 2 * pt * qt
  
  Hsn <-  2 * pn * qn
  Hs <- rowSums(Hsn) / ncol(Hsn)
  
  Fst <- (Ht - Hs) / Ht
  
  Fst[is.infinite(Fst)] <- 0
  Fst[is.nan(Fst)] <- 0
  
  summary(Fst)
  
  
  
  ## PLOT
  # Rectangles for 2R and 3R
  snp_table$Fst <- Fst
  
  rect_data <- data.frame(xmin = c(min(snp_table$ABS_POS[snp_table$CHROM == "2R"]),
                                   min(snp_table$ABS_POS[snp_table$CHROM == "3R"])),
                          xmax = c(max(snp_table$ABS_POS[snp_table$CHROM == "2R"]),
                                   max(snp_table$ABS_POS[snp_table$CHROM == "3R"])),
                          ymin = c(0,0),
                          ymax = c(0.7, 0.7),
                          col = c("grey", "grey"))
  
  axis_set <-
    c("2L" = (max(snp_table[snp_table$CHROM == "2L", ]$ABS_POS) + 
                min(snp_table[snp_table$CHROM == "2L", ]$ABS_POS))/2,
      "2R" = (max(snp_table[snp_table$CHROM == "2R", ]$ABS_POS) + 
                min(snp_table[snp_table$CHROM == "2R", ]$ABS_POS))/2,
      "3L" = (max(snp_table[snp_table$CHROM == "3L", ]$ABS_POS) + 
                min(snp_table[snp_table$CHROM == "3L", ]$ABS_POS))/2,
      "3R" = (max(snp_table[snp_table$CHROM == "3R", ]$ABS_POS) + 
                min(snp_table[snp_table$CHROM == "3R", ]$ABS_POS))/2,
      "X"  = (max(snp_table[snp_table$CHROM == "X",  ]$ABS_POS) + 
                min(snp_table[snp_table$CHROM == "X",  ]$ABS_POS))/2)
  
  p <- ggplot(snp_table, aes(x = ABS_POS, y = Fst)) +
    geom_line() +
    scale_x_continuous(labels = names(axis_set), breaks = axis_set) +
    geom_rect(data = rect_data, aes(xmin = xmin,
                                    xmax = xmax,
                                    ymin = ymin,
                                    ymax = ymax,
                                    fill = "grey"),
              inherit.aes = FALSE,
              alpha = 0.5) +
    scale_fill_identity() +
    ggtitle("Fst - O1-10 gen 20") +
    theme_bw() +
    labs(x = "Genomic position", y = "Fst")
  
  return(p)
}








## PLOT - window
# Rectangles for 2R and 3R
rect_data <- data.frame(xmin = c(min(snp_table$ABS_POS[snp_table$CHROM == "2R"]),
                                 min(snp_table$ABS_POS[snp_table$CHROM == "3R"])),
                        xmax = c(max(snp_table$ABS_POS[snp_table$CHROM == "2R"]),
                                 max(snp_table$ABS_POS[snp_table$CHROM == "3R"])),
                        ymin = c(0,0),
                        ymax = c(max(Fst)+0.1, max(Fst)+0.1),
                        col = c("grey", "grey"))

axis_set <-
  c("2L" = (max(snp_table[snp_table$CHROM == "2L", ]$ABS_POS) + 
              min(snp_table[snp_table$CHROM == "2L", ]$ABS_POS))/2,
    "2R" = (max(snp_table[snp_table$CHROM == "2R", ]$ABS_POS) + 
              min(snp_table[snp_table$CHROM == "2R", ]$ABS_POS))/2,
    "3L" = (max(snp_table[snp_table$CHROM == "3L", ]$ABS_POS) + 
              min(snp_table[snp_table$CHROM == "3L", ]$ABS_POS))/2,
    "3R" = (max(snp_table[snp_table$CHROM == "3R", ]$ABS_POS) + 
              min(snp_table[snp_table$CHROM == "3R", ]$ABS_POS))/2,
    "X"  = (max(snp_table[snp_table$CHROM == "X",  ]$ABS_POS) + 
              min(snp_table[snp_table$CHROM == "X",  ]$ABS_POS))/2)

ggplot(fst_window, aes(x = start_pos, y = mean_fst)) +
  scale_x_continuous(labels = names(axis_set), breaks = axis_set) +
  #facet_wrap(~ CHROM, scales = "free_x") +
  geom_rect(data = rect_data, aes(xmin = xmin,
                                  xmax = xmax,
                                  ymin = ymin,
                                  ymax = ymax,
                                  fill = "grey"),
            inherit.aes = FALSE,
            alpha = 0.5) +
  scale_fill_identity() +
  geom_line() +
  theme_bw() +
  labs(x = "Genomic position", y = "Fst")







## Correlation between Fst and -log10(p-value)
cmh_pvals <- readRDS("results/cmh_pvals.rds")

# Add the CMH vals to the snp_table so we can filter them all together
cmh_pvals$ABS_POS <- NULL
cmh_pvals$CHROM <- NULL
snp_table_cmh <- cbind(snp_table, cmh_pvals)

threshold <- 1e-100

column <- snp_table_cmh$cmh_adapted_o01_vs_o20_scaled

significant_snp_table <- snp_table_cmh[which(p.adjust(column,  method = "fdr") < threshold),]

summary(significant_snp_table$Fst)
summary(significant_snp_table$cmh_adapted_fdr_o01_vs_o20_scaled)

sig_snp_table_complete <- significant_snp_table[complete.cases(significant_snp_table),]

cor(sig_snp_table_complete$Fst, -log10(sig_snp_table_complete$cmh_adapted_fdr_o01_vs_o20_scaled))

plot(sig_snp_table_complete$Fst, -log10(sig_snp_table_complete$cmh_adapted_fdr_o01_vs_o20_scaled), xlab = expression("F"["st"]), ylab = expression("-log"[10]*"(p-values)"))
