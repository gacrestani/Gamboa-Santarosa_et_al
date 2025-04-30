library(ggplot2)     


df <- read.csv("~/Downloads/pca_data.csv")
df$Sample <- NULL
df$af <- NULL
df$Treatment_Rep <- paste0(df$Treatment, "_Rep", df$Replicate)
df$Timepoint <- as.factor(df$Timepoint)
df$Replicate <- as.character(df$Replicate)

treatment_colors <- c(
  "A" = "blue",
  "F" = "orange",
  "L" = "yellow",
  "H" = "red"
)

treatment_order <- c("A", "H", "F", "L")

# Plot
ggplot(df, aes(x = Timepoint, y = PC5, group = Treatment_Rep, color = Treatment)) +
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