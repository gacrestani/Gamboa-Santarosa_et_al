source("scripts/utils.R")


snp_table <- readRDS("data/processed/processed_snps_shahrestani.rds")
freq <- GetFreq(snp_table)

pca_df <- t(as.matrix(freq)) # Transpose the data frame to fit prcomp standard
pca <- prcomp(pca_df)

pca_data <- data.frame(sample=gsub("alt_|N_", "", colnames(freq)),
                       X = pca$x[,1],
                       Y = pca$x[,2],
                       Z = pca$x[,3],
                       row.names = NULL) # Select first three Principal Components

pca_data$Population <- factor(gsub("([A-Z]+)_rep.._(gen..)", "\\1_\\2", pca_data$sample))
pca_data$variance <- pca$sdev^2
pca_data$variance_percentage <- round(pca_data$variance / sum(pca_data$variance)*100, 2)

clustering_result <- kmeans(pca_df, centers = 3, nstart = 25)
pca_data$Cluster <- as.factor(clustering_result$cluster)

pca_data$Generation <- gsub(".*_(gen\\d+)", "\\1", pca_data$sample)
pca_data$Generation <- gsub("gen", "", pca_data$Generation)
pca_data$Generation <- gsub("01", "1", pca_data$Generation)
pca_data$Generation <- gsub("20|56", "20|56", pca_data$Generation)

pca_data$Treatment <- gsub("_gen..", "", pca_data$Population)

pca_plot <- ggplot(
  data = pca_data,
  aes(
    x=X,
    y=Y,
    color = Treatment,
    shape = Generation)) +
  geom_point(size=3, alpha=0.75) +
  xlab(paste0("PC1 - ", as.character(pca_data$variance_percentage[1]), "%")) +
  ylab(paste0("PC2 - ", as.character(pca_data$variance_percentage[2]), "%")) +
  scale_color_manual(values = c("#003F88", "#4F9DF2",      
                                "#8E0000", "#FF3B30")) +
  theme_bw() +
  stat_ellipse(aes(group = Cluster), level=0.95, linetype = 2)

# Panel 2
pca_data <- pca_data %>%
  mutate(Treatment_Rep = str_split_fixed(sample, "_", 3)[,1:2] %>%
           apply(1, paste, collapse = "_"))

treatment_order <- c("OBO", "OB", "nBO", "nB")
treatment_colors <- c("OBO" = "#003F88", "OB" = "#4F9DF2", "nBO" = "#8E0000", "nB" = "#FF3B30")

pc1 <- ggplot(pca_data, aes(x = Generation, y = X, group = Treatment_Rep, color = Treatment)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  scale_color_manual(
    values = treatment_colors,
    breaks = treatment_order,
    name = "Treatment"
  ) +
  theme_bw() +
  labs(
    #title = "PC1 Across Timepoints by Treatment",
    x = "Generations", y = paste0("PC1 - ", as.character(pca_data$variance_percentage[1]), "%")
  )  +
  theme(
    legend.position = "none"
  )


pca_combined <- ggarrange(
  pca_plot,
  pc1,
  ncol = 2,
  nrow = 1,
  labels = c("A", "B"),
  common.legend = TRUE,
  legend = "right"
)

ggsave(
  filename = "results/figure_3_pca.png",
  plot = pca_combined,
  width = 10,
  height = 5,
  dpi = 300
)
