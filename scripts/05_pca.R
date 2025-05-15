source("scripts/functions.R")

snp_table_shahrestani <- 
  readRDS("data/processed/processed_snps_abcd_shahrestani.rds")


df <- pca
#PlotPca(pca, label = FALSE)


library(ggplot2)     
library(stringr)

df <- PreparePca(snp_table_shahrestani)
#df <- read.csv("~/Downloads/pca_data.csv")
#df$Sample <- NULL
#df$af <- NULL

df <- df %>%
  mutate(Treatment_Rep = str_split_fixed(sample, "_", 3)[,1:2] %>%
           apply(1, paste, collapse = "_"))


df$Timepoint <- str_extract(df$Population, "\\d+")
df$Timepoint[df$Timepoint %in% c(20, 56)] <- 02

df$Treatment <- str_extract(df$Population, "^[^_]+")

treatment_colors <- c(
  "OBO" = "blue",
  "OB" = "#377EB8",
  "nBO" = "darkred",
  "nB" = "#E41A1C"
)

treatment_order <- c("OBO", "OB", "nBO", "nB")

# Plot
ggplot(df, aes(x = Timepoint, y = Y, group = Treatment_Rep, color = Treatment)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  scale_color_manual(
    values = treatment_colors,
    breaks = treatment_order,
    name = "Treatment"
  ) +
  theme_minimal() +
  labs(
    title = "PC1 Across Timepoints by Treatment",
    x = "Timepoint", y = "PC1"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

