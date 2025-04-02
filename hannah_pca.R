# Hannah's PCA code
# First, load all the functions you will need:
library(ggplot2)
library(dplyr)

FilterMinAndMaxCov <-
  function(
    snp_table,
    min_cov = 30,
    max_cov = 500) {
    
    snp_table <- snp_table %>%
      filter(if_all(starts_with("N_"), ~ . >= min_cov & . <= max_cov))
    
    return(snp_table)
  }

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
    # Hannah, change this to a code that gets your populations instead of mine (like F_01_12)
    # ChatGPT is your friend here!
    pca_data$Population <- factor(
      gsub("([A-Z]+)_rep.._(gen..)",
           "\\1_\\2",
           pca_data$sample))
    
    pca_data$variance <- pca$sdev^2 # Create var column
    pca_data$variance_percentage <-
      round(pca_data$variance / sum(pca_data$variance)*100, 2)
    
    clustering_result <- kmeans(pca_df, centers = 3, nstart = 25) # You can opt to do the clustering or not. If you decide not to, just delete this and the next line
    pca_data$Cluster <- as.factor(clustering_result$cluster) # This line also needs to be deleted
    
    return(pca_data)
  }

PlotPca <-
  function(
    pca_data,
    label = TRUE,
    title = NULL) {
    
    # If label = TRUE, all data points will be labeled independently
    # If label = FALSE, data points will be colored by population
    
    # Pick a palette that fits your data. You may want to change this.
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
                 shape = Cluster)) + # Delete shape = Cluster (and the final comma in the previous line) if not using k-means clustering
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

# Load your SNP table as an object called snp_table
# Your code may look like the following:
snp_table <- read.table("", header = TRUE) # Insert your table location in between ""
snp_table[snp_table == "."] <- 0 # For some reason, we have some dots in our table. Let's change them to 0s
snp_table[, 6:ncol(snp_table)] <- lapply(snp_table[, 6:ncol(snp_table)], as.numeric) # Make sure all alt_ and N_ columns are numeric
snp_table <- FilterMinAndMaxCov(snp_table, min_cov = 5, max_cov = 500) # This is a function I wrote to filter min and max coverage

# Now we will run the PCA and create the plot
pca <- PreparePca(snp_table)
PlotPca(pca, label = FALSE)


