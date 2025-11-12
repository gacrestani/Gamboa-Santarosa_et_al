# utils.R
# This script works like a library (in fact, it should become a library in the future).
# It contains all the functions I require to perform the analyses, as well as the libraries I need for each step.

# Libraries --------------------------------------------------------------------
library(parallel)
library(poolSeq)
library(ACER)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(ggrepel)
library(biomaRt)

options(scipen=999) # Disable scientific notation

'%!in%' <- function(x,y)!('%in%'(x,y))

FilterSamples <- 
  function(
    snp_table,
    treatment1,
    gen1,
    treatment2,
    gen2) {
  
  sample1 <- paste(treatment1, "_rep.._gen", gen1, sep = "")
  sample2 <- paste(treatment2, "_rep.._gen", gen2, sep = "")
  
  snp_table_filtered <- 
    snp_table[grep(paste(sample1, "|", sample2, sep = ""),
                   colnames(snp_table))]
  
  return(snp_table_filtered)
}

CountReplicates <-
  function(
    snp_table) {
  
  freq <- GetFreq(snp_table)
  replicates <- unique(gsub("[a-zA-Z]+_[a-zA-Z]+_rep(..)_gen..",
                            "\\1",
                            colnames(freq)))
  
  return(replicates)
}

GetFreq <-
  function(
    snp_table) {
  
  alt <- snp_table[,grep("^alt_", colnames(snp_table))]
  cov <- snp_table[,grep("^N_", colnames(snp_table))]
  
  freq <- alt / cov
  
  return(freq)
  }

GetCov <-
  function(
    snp_table) {
    
    cov <- snp_table[,grep("^N_", colnames(snp_table))]
    
    return(cov)
  }


# 4 Plotting functions =========================================================


PlotCoverage <-
  function(
    snp_table,
    treatment,
    gen,
    mode = "normal") {
  
  # Mode can be "normal", "scaled", or "standardized"
  
  # Scaling function
  scaleCol <- function(val) {
    scaledCol <- val / mean(val)
    return(scaledCol)
  }
  
  # Standardize function
  standardize = function(x){
    z <- (x - mean(x)) / sd(x)
    return( z)
  }
  
  cov <- snp_table[grep(paste("N_", treatment, "_rep.._gen", gen, sep=""), colnames(snp_table))]
  
  if (mode == "scaled") {
    cov2 <- as.data.frame(apply(cov, 2, scaleCol))
  }
  
  if (mode == "standardized") {
    cov3 <- as.data.frame(apply(cov, 2, standardize))
  }
  
  plots <- list()
  for(i in 1:ncol(cov)){
    plots[[i]] <- manhplot(snp_table, Y = cov[,i], title = paste(names(cov[i]), "_", mode, sep = ""))
  }
  do.call(grid.arrange, c(plots, ncol=3, nrow=2))
}

GetAllelicTrajectoryPlot <-
  function(
    snp_table,
    snp,
    treatments = c("OB", "OBO", "nB", "nBO"),
    cex = 2.5) {
  
  # Define samples and filter freq dataset to contain only the interesting SNP
  snp_chrom <- gsub(":", "", unlist(strsplit(snp, ":"))[1])
  snp_pos <- as.numeric(unlist(strsplit(snp, ":"))[2])
  
  #samples <- unique(gsub("alt_", "", gsub("N_", "", colnames(data[grep("alt_|N_", colnames(data),)]))))
  freq <- snp_table[,grep("^alt_", colnames(snp_table))]/snp_table[,grep("^N_", colnames(snp_table))]
  freq$POS <- snp_table$POS
  freq$CHROM <- snp_table$CHROM
  freq_filtered <- freq[snp_table$POS == snp_pos & snp_table$CHROM == snp_chrom, ]
  
  
  # Transpose the dataframe to make manipulation easier
  freq_filtered_t <- as.data.frame(t(freq_filtered))
  rownames(freq_filtered_t) <- gsub("alt_", "", rownames(freq_filtered_t))
  
  # Set the axis based on what treatments I am interested in
  if (any(c("CBO", "CB", "nBO", "nB") %in% treatments)) {axis_set = c("01" = 1, "56" = 20)}
  if (any(c("EBO", "EB", "OBO", "OB") %in% treatments)) {axis_set = c("01" = 1, "20" = 20)}
  if (any(c("CBO", "CB", "nBO", "nB") %in% treatments) & any(c("EBO", "EB", "OBO", "OB") %in% treatments)) {axis_set = c("01" = 1, "20/56" = 20)}
  
  # Set the title
  title = paste("SNP ", freq_filtered_t["CHROM", ], ":", format(as.numeric(freq_filtered_t["POS", ]), big.mark = ","), sep = "")
  
  colors <- c("EB" = "blue",
              "OB" = "blue",
              "EBO" = "darkblue",
              "OBO" = "darkblue",
              "CB" = "red",
              "nB" = "red",
              "CBO" = "darkred",
              "nBO" = "darkred")
  
  strokes <- c("EB" = 2,
               "OB" = 2,
               "EBO" = 1,
               "OBO" = 1,
               "CB" = 2,
               "nB" = 2,
               "CBO" = 1,
               "nBO" = 1)
  
  # Initialize a clean plot with defined lims
  plot(1,
       type = "n",
       xlab = "", 
       ylab = "",
       xlim = c(1, 20),  
       ylim = c(0, 1),
       xaxt = "n",
       main = title,
       cex.axis = cex
  ) 
  
  # Grabs the frequency for all desired treatments
  for (treatment in treatments) {
    plotting_df <- cbind(
      freq_filtered_t[grep(paste(treatment, "_rep0._gen01", sep=""), rownames(freq_filtered_t)), 1],
      freq_filtered_t[grep(paste(treatment, "_rep.._gen20|", treatment, "_rep.._gen56", sep=""), rownames(freq_filtered_t)), 1]
    )
    
    # Gets color and stroke for treatments
    color = colors[treatment]
    stroke = strokes[treatment]
    
    # Adds a line for each treatment
    for (line in 1:nrow(plotting_df)) {
      lines(c(1,20), plotting_df[line, ], type = "b", col=color, lty=stroke)
    }
    axis(1, at=axis_set, labels=names(axis_set), cex.axis = cex)
  }
  
  # Gets the legend. Commented for now as I don't require it
  # legend(x = "center",
  #        legend = treatments,
  #        lty = strokes[treatments],
  #        col = colors[treatments],
  #        cex = cex/2)
  
  plot <- recordPlot()
  
  return(plot)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 5 PCA ------------------------------------------------------------------------
PreparePca <- 
  function(
    snp_table) {
  
  freq <- GetFreq(snp_table)
  pca_df <- t(as.matrix(freq)) # Transpose the data frame to fit prcomp standard
  pca <- prcomp(pca_df)
  
  pca_data <- data.frame(sample=gsub("alt_|N_", "", colnames(freq)),
                         X = pca$x[,1],
                         Y = pca$x[,2],
                         Z = pca$x[,3]) # Select first three Principal Components
  
  # Create grouping variable
  pca_data$Population <- factor(
    gsub("([A-Z]+)_rep.._(gen..)",
         "\\1_\\2",
         pca_data$sample))
  
  pca_data$variance <- pca$sdev^2 # Create var column
  pca_data$variance_percentage <-
    round(pca_data$variance / sum(pca_data$variance)*100, 2)
  
  clustering_result <- kmeans(pca_df, centers = 3, nstart = 25)
  pca_data$Cluster <- as.factor(clustering_result$cluster)
  
  return(pca_data)
}

PlotPca <-
  function(
    pca_data,
    label = TRUE,
    title = NULL) {
  
  # If label = TRUE, all data points will be labeled independently
  # If label = FALSE, data points will be colored by population
  
  # Pick a palette that fits your data
  rbpalette=c("red",
              "darkred",
              "magenta",
              "darkorchid4",
              "blue",
              "darkblue",
              "green",
              "darkgreen")
  
  # This one has the text label for all the samples
  pca_plot <-
    ggplot(data = pca_data,
           aes(x=X,
               y=Y,
               label=sample,
               color=Population,
               shape = Cluster)) +
    geom_point(size = 3) +
    scale_color_manual(values = rbpalette) +
    xlab(paste("PC1 - ", pca_data$variance_percentage[1], "%", sep="")) +
    ylab(paste("PC2 - ", pca_data$variance_percentage[2], "%", sep="")) +
    ggtitle(title)  +
    {if (label) geom_text_repel()} +
    {if (label) theme(legend.position="none")} +
    theme_bw()
  
  return(pca_plot)
  }
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 6 GO Analysis ----------------------------------------------------------------
FilterOutFixedSnps <-
  function(
    snp_table) {
  
  o_type_samples <- "alt_OB_rep0._gen01|alt_OBO_rep0._gen01"
  
  selected_columns <-
    colnames(snp_table)[grep(o_type_samples, colnames(snp_table))]
  
  rows_to_exclude <-
    apply(freq[selected_columns],
          1,
          function(row) {any(row == 0) || any(row == 1)})
  
  filtered_snp_table <- snp_table[!rows_to_exclude, ]
  
  rows_lost <- sum(rows_to_exclude)
  
  message <-
    paste("Number of SNPs discarded:", format(rows_lost, big.mark = ","))
  
  print(message)
  
  return(filtered_snp_table)
}

DiagnoseSnps <-
  function(
    filtered_snp_table) {
  for (i in 1:nrow(filtered_snp_table)) {
    
    GetAllelicTrajectoryPlot(snp_table = filtered_snp_table,
                             snp = paste(filtered_snp_table$CHROM[i], ":", filtered_snp_table$POS[i], sep = ""),
                             treatments = c("OBO", "OB"))
    
    cat("Row ", i, "\n")
    
    # Prompt user input
    response <- readline(prompt = "Press enter to continue. Type anything to exit.\n")
    
    # Check response
    if (toupper(response) != "") {
      cat("Exiting.\n")
      break
    }
    
    cat("Continuing to the next SNP.\n")
  }
}

GetGeneList <-
  function(
    GO_dataframe,
    peaks,
    ensembl,
    window = 10000) {
    
    gene_list <- c()
    
    for (i in 1:length(peaks)) {
      chromosome <- GO_dataframe$CHROM[GO_dataframe$coordinate == peaks[i]]
      pos <- GO_dataframe$POS[GO_dataframe$coordinate == peaks[i]]
      
      result <- getBM(
        attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
        filters = c("chromosome_name", "start", "end"),
        values = list(chromosome, pos-window, pos+window),
        mart = ensembl
      )
      
      print(result)
      gene_list <- c(gene_list, result$external_gene_name)
    }
    
    return(gene_list)
    
  }
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~