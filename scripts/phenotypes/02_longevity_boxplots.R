library(readxl)
library(tidyverse)

input_df <- read_excel("data/phenotypic_data/raw/Longevity Gen 20.xlsx") %>%
  select(-sign, -Notes)

days <- as.numeric(colnames(input_df)[3:ncol(input_df)])
days <- c(0, days[-1] - days[1])
d_days <- paste0("d", days)
colnames(input_df) <- c("Cage", "Sex", d_days)

longevity_df <- input_df[d_days]
pop_sizes <- rowSums(longevity_df)
longevity_df <- sweep(longevity_df, MARGIN = 2, days, `*`)

mean_cage_longevity <- rowSums(longevity_df) / pop_sizes
summary_df <- cbind(input_df[, 1:2], mean_cage_longevity)

summary_df <- summary_df %>%
  mutate(Population = str_remove(Cage, " \\(.*\\)")) %>%
  group_by(Population, Sex) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

calc_df <- summary_df %>%
  mutate(
    mean_pop_longevity = mean_cage_longevity,
    Population = summary_df$Population,
    Regime = case_when(
      grepl("EB|EBO", Population) ~ "O-type",
      grepl("CB|CBO", Population) ~ "B-type",
      TRUE ~ NA_character_
    ),
    Ancestry = case_when(
      grepl("CBO|EBO", Population) ~ "BO",
      grepl("CB|EB", Population) ~ "B",
      TRUE ~ NA_character_
    ),
    Sex = case_when(
      grepl("F", Sex) ~ "Female",
      grepl("M", Sex) ~ "Male"
    ),
    Replicate = as.factor(gsub("\\D", "", Population)),
    Treatment = case_when(
      grepl("EBO", Population) ~ "OBO",
      grepl("EB", Population) ~ "OB",
      grepl("CBO", Population) ~ "nBO",
      grepl("CB", Population) ~ "nB",
      TRUE ~ NA_character_
    )
  )

createPlot <- function(ancestry = FALSE) {
  p <- ggplot(calc_df, aes(x = Sex, y = mean_pop_longevity, fill = Regime)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
    theme_bw() +
    labs(
      #title = "Longevity",
      x = NULL,
      y = "Longevity (days)"
    ) +
    theme(legend.position = "none") 
  
  if (ancestry) {
    p <- p + facet_grid(~Ancestry)
  }
  
  p <- p + stat_compare_means(aes(group = Regime), label = "p.signif", method = "t.test")
  
  return(p)
}

longevity_boxplots <- createPlot()
longevity_boxplots_anc <- createPlot(ancestry = TRUE)

longevity_boxplots 
longevity_boxplots_anc
# Statistical analysis
# Difference in population mean longevity
lm_fit_longevity <- lm(mean_pop_longevity ~ Regime + Ancestry + Sex, data = calc_df)
summary(lm_fit_longevity)