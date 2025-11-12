# cmh_tests.R
# Runs the CMH tests on the SNP tables.
#
# inputs: snp_table.rds (x4)
# outputs: cmh_pvals.rds - a dataframe with the p-values for each test.

source("scripts/utils.R")

# Create used functions
GetNe <-
  function(snp_table,
           treatment1,
           gen1,
           treatment2,
           gen2,
           t = 20) {
    # Ne stands for effective population size. This function uses the methods of ACER to estimate Ne
    # Parameters of the estimateNe function are adjusted for my samples
    alt_cmh1 <- snp_table[grep(paste0("alt_", treatment1, "_rep.._gen", gen1),
                               colnames(snp_table))]
    alt_cmh2 <- snp_table[grep(paste0("alt_", treatment2, "_rep.._gen", gen2),
                               colnames(snp_table))]
    cov_cmh1 <- snp_table[grep(paste0("N_", treatment1, "_rep.._gen", gen1),
                               colnames(snp_table))]
    cov_cmh2 <- snp_table[grep(paste0("N_", treatment2, "_rep.._gen", gen2),
                               colnames(snp_table))]
    
    freq_cmh1 <- alt_cmh1 / cov_cmh1
    freq_cmh2 <- alt_cmh2 / cov_cmh2
    
    # Creates list of Ne and calculate values for each replicate
    Ne <- c()
    for (i in 1:ncol(freq_cmh1)) {
      estimated_Ne <- estimateNe(
        p0 = freq_cmh1[, i],
        pt = freq_cmh2[, i],
        cov0 = cov_cmh1[, i],
        covt = cov_cmh2[, i],
        t = t,
        method = c("P.planII"),
        poolSize = c(100, 100)
      )
      
      # In case I want to save these
      pop_name <- gsub(".*_([A-Z]+_rep\\d+).*", "\\1", colnames(freq_cmh1[i]))
      cat(pop_name, "estimated Ne:", estimated_Ne, "\n")
      Ne <- c(Ne, as.integer(unname(estimated_Ne)))
    }
    
    return(Ne)
  }

AdaptedCmhTest <-
  function(snp_table,
           treatment1,
           gen1,
           treatment2,
           gen2,
           t,
           Ne) {
    snp_table_filtered <-
      FilterSamples(
        snp_table = snp_table,
        treatment1 = treatment1,
        gen1 = gen1,
        treatment2 = treatment2,
        gen2 = gen2
      )
    
    replicates <- CountReplicates(snp_table = snp_table_filtered)
    
    freq <- GetFreq(snp_table = snp_table_filtered)
    
    cov <- snp_table_filtered[, grep("^N_", colnames(snp_table_filtered))]
    
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
        t = 20
      )
    
    pvals <-
      adapted.cmh.test(
        freq = as.matrix(freq),
        coverage = as.matrix(cov),
        Ne = Ne,
        gen = as.numeric(c(gen1, gen2)),
        repl = 1:length(replicates),
        poolSize = rep(c(100, 100), length(replicates)),
        mincov = 1,
        MeanStart = TRUE,
        IntGen = FALSE,
        TA = FALSE,
        order = 0,
        correct = FALSE,
        RetVal = 0
      )
    
    return(pvals)
  }

ClassicalCmhTest <- function(snp_table,
                             treatment1,
                             gen1,
                             treatment2,
                             gen2) {
  # Filter table
  snp_table_filtered <- FilterSamples(
    snp_table = snp_table,
    treatment1 = treatment1,
    gen1 = gen1,
    treatment2 = treatment2,
    gen2 = gen2
  )
  
  replicates <- CountReplicates(snp_table)
  n_reps <- length(replicates)
  
  # Precompute column names
  alt1_cols <- paste0("alt_", treatment1, "_rep", replicates, "_gen", gen1)
  alt2_cols <- paste0("alt_", treatment2, "_rep", replicates, "_gen", gen2)
  n1_cols   <- paste0("N_", treatment1, "_rep", replicates, "_gen", gen1)
  n2_cols   <- paste0("N_", treatment2, "_rep", replicates, "_gen", gen2)
  
  # Pre-subset numeric matrices once
  alt1 <- as.matrix(snp_table_filtered[, alt1_cols])
  alt2 <- as.matrix(snp_table_filtered[, alt2_cols])
  n1   <- as.matrix(snp_table_filtered[, n1_cols])
  n2   <- as.matrix(snp_table_filtered[, n2_cols])
  
  # Run CMH test per SNP (row)
  num_cores <- max(1, parallel::detectCores(logical = TRUE) - 1)
  
  p_values <- parallel::mclapply(seq_len(nrow(snp_table_filtered)), function(i) {
    vals <- rbind(alt1[i, ], alt2[i, ], n1[i, ], n2[i, ])
    
    if (anyNA(vals))
      return(1)  # immediately skip
    
    # Build 2x2xreplicate array
    a <- array(c(vals[1, ], vals[3, ] - vals[1, ], vals[2, ], vals[4, ] - vals[2, ]),
               dim = c(2, 2, n_reps))
    
    # Try running the test
    out <- tryCatch(
      stats::mantelhaen.test(a)$p.value,
      error = function(e)
        1
    )
    if (is.nan(out))
      out <- 1
    out
  }, mc.cores = num_cores)
  
  p <- unlist(p_values)
  return(p)
  
}

# Read snp tables
snp_table_shahrestani_scaled <-
  readRDS("data/processed/processed_snps_shahrestani_scaled.rds")

snp_table_regimes_scaled <-
  readRDS("data/processed/processed_snps_regimes_scaled.rds")

# O
adapted_o20_vs_o01 <- AdaptedCmhTest(
  snp_table = snp_table_regimes_scaled,
  treatment1 = "O",
  gen1 = "01",
  treatment2 = "O",
  gen2 = "20",
  t = 20
)

adapted_obo20_vs_obo01 <- AdaptedCmhTest(
  snp_table = snp_table_shahrestani_scaled,
  treatment1 = "OBO",
  gen1 = "01",
  treatment2 = "OBO",
  gen2 = "20",
  t = 20
)

adapted_ob20_vs_ob01 <- AdaptedCmhTest(
  snp_table = snp_table_shahrestani_scaled,
  treatment1 = "OB",
  gen1 = "01",
  treatment2 = "OB",
  gen2 = "20",
  t = 20
)


# B
adapted_b56_vs_b01 <- AdaptedCmhTest(
  snp_table = snp_table_regimes_scaled,
  treatment1 = "B",
  gen1 = "01",
  treatment2 = "B",
  gen2 = "56",
  t = 56
)

adapted_nbo56_vs_nbo01 <- AdaptedCmhTest(
  snp_table = snp_table_shahrestani_scaled,
  treatment1 = "nBO",
  gen1 = "01",
  treatment2 = "nBO",
  gen2 = "56",
  t = 56
)

adapted_nb56_vs_nb01 <- AdaptedCmhTest(
  snp_table = snp_table_shahrestani_scaled,
  treatment1 = "nB",
  gen1 = "01",
  treatment2 = "nB",
  gen2 = "56",
  t = 56
)


# O vs B late generation
classical_o20_vs_b56 <- ClassicalCmhTest(
  snp_table = snp_table_regimes_scaled,
  treatment1 = "O",
  gen1 = "20",
  treatment2 = "B",
  gen2 = "56"
)

classical_obo20_vs_nbo56 <- ClassicalCmhTest(
  snp_table = snp_table_shahrestani_scaled,
  treatment1 = "OBO",
  gen1 = "20",
  treatment2 = "nBO",
  gen2 = "56"
)

classical_ob20_vs_nb56 <- ClassicalCmhTest(
  snp_table = snp_table_shahrestani_scaled,
  treatment1 = "OB",
  gen1 = "20",
  treatment2 = "nB",
  gen2 = "56"
)



# Adapted.cmh.test can't be performed on O vs B because it is not time-series data

cmh_pvals <- as.data.frame(cbind(
  adapted_o20_vs_o01     = adapted_o20_vs_o01,
  adapted_obo20_vs_obo01 = adapted_obo20_vs_obo01,
  adapted_ob20_vs_ob01   = adapted_ob20_vs_ob01,
  adapted_b56_vs_b01     = adapted_b56_vs_b01,
  adapted_nbo56_vs_nbo01 = adapted_nbo56_vs_nbo01,
  adapted_nb56_vs_nb01   = adapted_nb56_vs_nb01,
  classical_o20_vs_b56     = classical_o20_vs_b56,
  classical_obo20_vs_nbo56 = classical_obo20_vs_nbo56,
  classical_ob20_vs_nb56   = classical_ob20_vs_nb56
))


# FDR correct everything
fdr_cmh <- as.data.frame(lapply(cmh_pvals, function(p) p.adjust(p, method = "fdr")))
names(fdr_cmh) <- paste0(names(fdr_cmh), "_fdr")

all_cmh <- cbind(
  snp_table_shahrestani_scaled[c("CHROM", "POS", "REF", "ALT", "ABS_POS")],
  cmh_pvals,
  fdr_cmh)

# Save results
saveRDS(all_cmh, "results/all_cmh.rds")