
# functions.R
# This script works like a library (in fact, it should become a library in the future).
# It contains all the functions I require to perform the analyses, as well as the libraries I need for each step.
#
# inputs: none
# outputs: functions are loaded into the environment

#install.packages("devtools")
#devtools::install_github("ThomasTaus/poolSeq")
#devtools::install_github("MartaPelizzola/ACER")
#BiocManager::install("biomaRt")

# Libraries --------------------------------------------------------------------
library(parallel)
library(poolSeq)
library(ACER)
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(ggrepel)
library(biomaRt)

options(scipen=999) # Disable scientific notation

'%!in%' <- function(x,y)!('%in%'(x,y))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 0 Data Manipulation functions ------------------------------------------------
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

# 1 Data Preparation ===========================================================
ReadSnpTable <-
  function(
    path = "data/snp_tables/filtered_snps_abcd.txt",
    mode = "shahrestani") {
  
  snp_table <- data.table::fread(file = path, header = TRUE)
  snp_table <- as.data.frame(snp_table)
  
  # Replaces "." with 0 and converts columns 6 to the end (i.e. those with data) to numeric.
  # TO DO: why are we getting those dots?
  snp_table[snp_table == "."] <- 0
  snp_table[, 6:ncol(snp_table)] <- lapply(snp_table[, 6:ncol(snp_table)], as.numeric)
  
  # Determines the sample naming convention used in my project
  # default:     EBO, EB, CBO, CB
  # shahrestani: OBO, OB, nBO, nB
  # regimes:     O_1~10, B_1~10
  
  if (mode == "regimes") {
    colnames <- colnames(snp_table)
    colnames <- gsub("CBO", "B", colnames(snp_table))
    colnames <- gsub("CB_rep01", "B_rep06", colnames)
    colnames <- gsub("CB_rep02", "B_rep07", colnames)
    colnames <- gsub("CB_rep03", "B_rep08", colnames)
    colnames <- gsub("CB_rep04", "B_rep09", colnames)
    colnames <- gsub("CB_rep05", "B_rep10", colnames)
    
    colnames <- gsub("EBO", "O", colnames)
    colnames <- gsub("EB_rep01", "O_rep06", colnames)
    colnames <- gsub("EB_rep02", "O_rep07", colnames)
    colnames <- gsub("EB_rep03", "O_rep08", colnames)
    colnames <- gsub("EB_rep04", "O_rep09", colnames)
    colnames <- gsub("EB_rep05", "O_rep10", colnames)
    
    colnames(snp_table) <- colnames
    
  } else if (mode == "shahrestani") {
    colnames <- colnames(snp_table)
    
    colnames <- gsub("EB", "OB", colnames)
    colnames <- gsub("CB", "nB", colnames)
    
    colnames(snp_table) <- colnames
  }
  
  return(snp_table)
}

AddAbsPosToSnpTable <-
  function(
    snp_table) {
  
  # Drosophila chromosome size data from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4352887/ (table 1)
  chromosome_sizes <- c("2L" = 23513712,
                        "2R" = 25286936,
                        "3L" = 28110227,
                        "3R" = 32079331,
                        "X" = 23542271) 
  
  snp_table <- snp_table %>%
    mutate(ABS_POS = case_when(CHROM == "2L" ~ POS,
                               CHROM == "2R" ~ POS + sum(chromosome_sizes[1]),
                               CHROM == "3L" ~ POS + sum(chromosome_sizes[1:2]),
                               CHROM == "3R" ~ POS + sum(chromosome_sizes[1:3]),
                               CHROM == "X" ~ POS  + sum(chromosome_sizes[1:4])))
  
  return(snp_table)
}

FilterMAF <-
  function(
    snp_table,
    limit = 0.001) {
  
  # MAF stands for minor allele frequency
  # If all samples are below a MAF cutoff, the SNP is removed
  samples <- gsub("N_", "", colnames(snp_table)[grepl("^N_", colnames(snp_table))])
  alt_sums <- rowSums(snp_table[paste("alt_", samples, sep="")])
  total_sums <- rowSums(snp_table[paste("N_", samples, sep="")])
  freq_rowsums <- alt_sums / total_sums
  snp_table <- snp_table[freq_rowsums >= limit & freq_rowsums <= (1 - limit),]
  
  return(snp_table)
}

FilterMinAndMaxCov <-
  function(
    snp_table,
    min_cov = 30,
    max_cov = 500) {
  
  snp_table <- snp_table %>%
    filter(if_all(starts_with("N_"), ~ . >= min_cov & . <= max_cov))
  
  return(snp_table)
}

FilterChromosomes <-
  function(
    snp_table = snp_table,
    chromosomes = c("2L", "2R", "3L", "3R", "X")) {
  
  snp_table <- snp_table[snp_table$CHROM %in% chromosomes, ]
  
  return(snp_table)
}

ReadAndPrepare <-
  function(
    path = "data/snp_tables/filtered_snps_abcd.txt",
    mode = "shahrestani",
    limit = 0.001,
    min_cov = 30,
    max_cov = 500,
    chromosomes = c("2L", "2R", "3L", "3R", "X")) {
  
  snp_table <- ReadSnpTable(path, mode)
  snp_table <- AddAbsPosToSnpTable(snp_table)
  snp_table <- FilterMAF(snp_table, limit)
  snp_table <- FilterMinAndMaxCov(snp_table, min_cov, max_cov)
  snp_table <- FilterChromosomes(snp_table, chromosomes)
  
  return(snp_table)
  
}

ScaleSnpTable <-
  function(snp_table) {
  
  # Identify columns
  alt_cols <- grep("^alt_", names(snp_table), value = TRUE)
  cov_cols <- grep("^N_", names(snp_table), value = TRUE)
  
  # Scale alt_ columns
  for (i in seq_along(alt_cols)) {
    snp_table[[alt_cols[i]]] <- snp_table[[alt_cols[i]]] * (100 / snp_table[[cov_cols[i]]])
  }
  
  snp_table[cov_cols] <- 100
  
  return(snp_table)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 2 CMH Test -------------------------------------------------------------------
ClassicalCmhTest <-
  function(
    snp_table,
    treatment1,
    gen1 = "01",
    treatment2 = treatment1,
    gen2 = "20") {
  
  # Force evaluation of treatment2
  force(treatment2)
  
  # Changes gen2 if dealing with nBO or nB or B
  if (treatment2 == "nBO" | treatment2 == "nB" | treatment2 == "B") {
    gen2 <- "56"
  }
  
  snp_table_filtered <-
    FilterSamples(
      snp_table = snp_table,
      treatment1 = treatment1,
      gen1 = gen1,
      treatment2 = treatment2,
      gen2 = gen2)
  
  replicates <- CountReplicates(snp_table)
  
  # Picks the columns we will use to build the matrix
  alt1_cols <- paste0("alt_", treatment1, "_rep", replicates, "_gen", gen1)
  alt2_cols <- paste0("alt_", treatment2, "_rep", replicates, "_gen", gen2)
  n1_cols <- paste0("N_", treatment1, "_rep", replicates, "_gen", gen1)
  n2_cols <- paste0("N_", treatment2, "_rep", replicates, "_gen", gen2)
  
  # Coerces the dataset into a matrix of 2x2xreplicates dimensions
  # This uses the parallel package to run the loops in parallel
  p_list <- mclapply(1:nrow(snp_table_filtered), function(line) {
    
    vals <- rbind(
      as.numeric(snp_table_filtered[line, alt1_cols]),
      as.numeric(snp_table_filtered[line, alt2_cols]),
      as.numeric(snp_table_filtered[line, n1_cols]),
      as.numeric(snp_table_filtered[line, n2_cols])
    )
    
    # Check for missing values across all replicates
    if (any(is.na(vals))) return(NA_real_)
    
    # Reshape into a 3D array.
    # Each column of 'vals' corresponds to one replicate:
    # The 2x2 matrix is formed column-wise:
    #   - First column: (v1, v2)
    #   - Second column: (v3, v4)
    matrices <- array(vals, dim = c(2, 2, length(alt1_cols)))
    
    # Run the Mantel-Haenszel test if the matrix is valid
    test <- tryCatch(mantelhaen.test(matrices), error = function(e) NA_real_)
    test$p.value
  }, mc.cores = max(1, detectCores() - 1))
  
  # Convert results to a numeric vector
  p_list <- unlist(p_list)
  
  # Sometimes we get some NA p-values, which are converted to 1
  # This happens when all samples of one side of the comparison (i.e. all O_gen01 samples) have a frequency of 0
  p_list[is.nan(p_list)] <- 1
  
  return(p_list)
}

GetNe <- 
  function(
    snp_table,
    treatment1,
    gen1,
    treatment2,
    gen2,
    t = 20) {
  
  # Ne stands for effective population size. This function uses the methods of ACER to estimate Ne
  # Parameters of the estimateNe function are adjusted for my samples
  alt_cmh1 <- snp_table[grep(paste0("alt_", treatment1, "_rep.._gen", gen1), colnames(snp_table))]
  alt_cmh2 <- snp_table[grep(paste0("alt_", treatment2, "_rep.._gen", gen2), colnames(snp_table))]
  cov_cmh1 <- snp_table[grep(paste0("N_", treatment1, "_rep.._gen", gen1), colnames(snp_table))]
  cov_cmh2 <- snp_table[grep(paste0("N_", treatment2, "_rep.._gen", gen2), colnames(snp_table))]
  
  freq_cmh1 <- alt_cmh1/cov_cmh1
  freq_cmh2 <- alt_cmh2/cov_cmh2
  
  # Creates list of Ne and calculate values for each replicate
  Ne <- c()
  for (i in 1:ncol(freq_cmh1)){
    estimated_Ne <- estimateNe(
      p0 = freq_cmh1[,i],
      pt = freq_cmh2[,i],
      cov0 = cov_cmh1[,i],
      covt = cov_cmh2[,i],
      t = t,
      method = c("P.planII"),
      poolSize = c(100, 100))
    
    # In case I want to save these
    pop_name <- gsub(".*_([A-Z]+_rep\\d+).*", "\\1", colnames(freq_cmh1[i]))
    cat(pop_name, "estimated Ne:", estimated_Ne, "\n")
    Ne <- c(Ne, as.integer(unname(estimated_Ne)))
  }
  
  return(Ne)
}

AdaptedCmhTest <-
  function(
    snp_table,
    treatment1,
    gen1 = "01",
    treatment2 = treatment1,
    gen2 = "20",
    t = 20,
    Ne) {
  
  # Force evaluation of treatment2
  force(treatment2)
  
  # Changes gen2 if dealing with nBO or nB or B
  if (treatment2 == "nBO" | treatment2 == "nB" | treatment2 == "B") {
    gen2 <- "56"
  }
  
  snp_table_filtered <-
    FilterSamples(
      snp_table = snp_table,
      treatment1 = treatment1,
      gen1 = gen1,
      treatment2 = treatment2,
      gen2 = gen2)
  
  replicates <- CountReplicates(snp_table = snp_table_filtered)
  
  freq <- GetFreq(snp_table = snp_table_filtered)
  
  cov <- snp_table_filtered[,grep("^N_", colnames(snp_table_filtered))]
  
  # Checks if both generations are the same. If so, differentiate them
  if (gen1 == gen2) {
    gen1 <- 01
    gen2 <- 02
  }
  
  # Calculates Ne
  Ne <-
    GetNe(
      snp_table = snp_table_filtered,
      treatment1 = treatment1,
      gen1 = gen1,
      treatment2 = treatment2,
      gen2 = gen2,
      t = t)
  
  pvals <- 
    adapted.cmh.test(
      freq = as.matrix(freq),
      coverage = as.matrix(cov),
      Ne = Ne,
      gen = as.numeric(c(gen1,gen2)),
      repl = 1:length(replicates),
      poolSize = rep(c(100, 100), length(replicates)),
      mincov = 1,
      MeanStart = TRUE,
      IntGen = FALSE,
      TA = FALSE,
      order = 0,
      correct = FALSE,
      RetVal = 0)
  
  return(pvals)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 3 Permutation Test ===========================================================
RunPermutationTestIteration <- 
  function(
    method = c("classic", "adapted"),
    filename,
    iter,
    snp_table,
    treatment1,
    gen1,
    treatment2,
    gen2,
    Ne = NULL) {
    
    # WARNING!
    # I think you don't have to run the permutation test if using the adapted 
    # method. Double check literature!
    
    # Assures a correct method is selected
    method <- match.arg(method)
    
    # Measure time needed to run one iteration
    start_time <- Sys.time()
    
    # Filter df to contain only interesting columns
    snp_table_filtered <-
      FilterSamples(snp_table,
                    treatment1 = treatment1,
                    gen1 = gen1,
                    treatment2 = treatment2,
                    gen2 = gen2)
    
    # Define possible choices for the sampling process
    choices_pattern_a <- paste("alt_", treatment1, "_rep.._gen", gen1, sep="")
    choices_pattern_b <- paste("alt_", treatment2, "_rep.._gen", gen2, sep="")
    
    choices_a <- 
      colnames(
        snp_table_filtered[,grep(choices_pattern_a, colnames(snp_table_filtered))]
      )
    
    choices_b <- colnames(
      snp_table_filtered[,grep(choices_pattern_b, colnames(snp_table_filtered))]
    )
    
    choices <- c(choices_a, choices_b)
    
    # Sample (randomize) choices
    choices_rand_alt <- sample(choices, length(choices))
    choices_rand_N <- gsub("alt_", "N_", choices_rand_alt)
    
    choices_rand <- c(rbind(choices_rand_alt, choices_rand_N))
    
    #colnames(snp_table_filtered) <- choices_rand
    snp_table_filtered <- snp_table_filtered[,choices_rand]
    
    if (method == "classic") {
      
      pval_list <- 
        ClassicalCmhTest(
          snp_table = snp_table,
          treatment1 = treatment1,
          gen1 = gen1,
          treatment2 = treatment2,
          gen2 = gen2)
      
    } else if (method == "adapted") {
      
      pval_list <- 
        AdaptedCmhTest(
          snp_table = snp_table_filtered,
          treatment1 = treatment1,
          gen1 = gen1,
          treatment2 = treatment2,
          gen2 = gen2,
          Ne = Ne)
      
    }
    
    pval_list[pval_list == 0] <- 1
    pval_list_log <- -log10(pval_list)
    max_pval <- max(pval_list_log)
    
    end_time <- Sys.time()
    end_time-start_time
    print(paste("Iter", iter, "- Time elapsed:", end_time-start_time, "seconds.", "p-val:", max_pval))
    
    fwrite(x = as.list(max_pval),
           file = filename,
           append=TRUE)
    
    return(max_pval)
  }

# 4 Plotting functions =========================================================
GetManhattanPlot <-
  function(
    my_dataframe,
    Y,
    permutation_pvals = NULL,
    percentage_significance = NULL,
    title = "Manhattan Plot",
    x_label = TRUE,
    y_label = "-log10(p)",
    palette = "blue",
    y_limit_up = 100,
    y_limit_down = 0) {
    
    # This huge function plots a good manhattan plot considering all the settings I have on the arguments.
    # Some of them have been tweaked over the time to reach what I think is an optimal graph for my use case
    axis_set <-
      c("2L" = (max(my_dataframe[my_dataframe$CHROM == "2L", ]$ABS_POS) + 
                  min(my_dataframe[my_dataframe$CHROM == "2L", ]$ABS_POS))/2,
        "2R" = (max(my_dataframe[my_dataframe$CHROM == "2R", ]$ABS_POS) + 
                  min(my_dataframe[my_dataframe$CHROM == "2R", ]$ABS_POS))/2,
        "3L" = (max(my_dataframe[my_dataframe$CHROM == "3L", ]$ABS_POS) + 
                  min(my_dataframe[my_dataframe$CHROM == "3L", ]$ABS_POS))/2,
        "3R" = (max(my_dataframe[my_dataframe$CHROM == "3R", ]$ABS_POS) + 
                  min(my_dataframe[my_dataframe$CHROM == "3R", ]$ABS_POS))/2,
        "X"  = (max(my_dataframe[my_dataframe$CHROM == "X",  ]$ABS_POS) + 
                  min(my_dataframe[my_dataframe$CHROM == "X",  ]$ABS_POS))/2)
    
    p <- ggplot(my_dataframe,
                aes(x = ABS_POS, 
                    y = Y,
                    color = as.factor(CHROM))) +
      geom_point(alpha = 0.75) +
      scale_x_continuous(labels = names(axis_set), breaks = axis_set) +
      {if (!x_label) scale_x_continuous(labels = NULL, breaks = axis_set)} +
      ylim(y_limit_down, y_limit_up) +
      {if (palette == "blue")
        scale_color_manual(values = rep(c("steelblue1", "steelblue4"), unique(length(axis_set))))} +
      {if (palette == "red")
        scale_color_manual(values = rep(c("firebrick1", "firebrick4"), unique(length(axis_set))))} +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = NULL, y = y_label) +
      theme_minimal() +
      theme(legend.position="none") +
      ggtitle(title)
    
    if (!is.null(permutation_pvals)) {
      qt90 <- quantile(as.numeric(permutation_pvals), probs = 0.90)
      qt95 <- quantile(as.numeric(permutation_pvals), probs = 0.95)
      qt99 <- quantile(as.numeric(permutation_pvals), probs = 0.99)
      qt999 <- quantile(as.numeric(permutation_pvals), probs = 0.999)
      
      p <- p +
        # 90% line
        geom_hline(yintercept = qt90, col = "darkgreen") +
        annotate("text", y=(qt90-2), x=max(my_dataframe$ABS_POS), label = "10%", col = "darkgreen") +
        # 95% line
        geom_hline(yintercept = qt95, col = "darkorange") +
        annotate("text", y=(qt95-2), x=max(my_dataframe$ABS_POS), label = "5%", col = "darkorange") +
        # 99% line
        geom_hline(yintercept = qt99, col = "burlywood4") +
        annotate("text", y=(qt99-2), x=max(my_dataframe$ABS_POS), label = "1%", col = "burlywood4") +
        # 99.9% line
        geom_hline(yintercept = qt999, col = "darkorchid") +
        annotate("text", y=(qt999-2), x=max(my_dataframe$ABS_POS), label = "0.1%", col = "darkorchid")
    }
    
    if (percentage_significance) {
      threshold <- 100 # or -log10(1e-100)
      p <- p +
        geom_hline(yintercept = threshold, col = "red")
      #annotate("text", y=(threshold-5), x=max(my_dataframe$ABS_POS), label = paste(percentage_significance*100, "% top SNPs"), col = "red")
    }
    
    return(p)
  }

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