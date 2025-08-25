library(readxl)
library(tidyverse)

# Use this to generate the same plot but with the intermediate generation data
input_df <- read_excel("data/phenotypic_data/raw/Longevity Gen 12 (cleaned).xlsx") %>%
  select(-Notes, - `Cage #`)

days <- as.numeric(colnames(input_df)[3:ncol(input_df)])
days <- c(0, days[-1] - days[1])
colnames(input_df) <- c("Cage", "Sex", paste0("d", days))

summary_df <- input_df %>%
  group_by(Cage, Sex) %>%
  summarise(across(where(is.numeric), ~ mean(replace(.x, is.na(.x), 0))), .groups = "drop")

count_matrix <- summary_df[3:ncol(summary_df)]
total_flies <- rowSums(count_matrix, na.rm = TRUE)

longevity_matrix <- sweep(count_matrix, MARGIN = 2, days, `*`)
total_longevity <- rowSums(longevity_matrix)

calc_df <- count_matrix %>%
  mutate(
    mean_longevity = total_longevity / total_flies,
    Cage = summary_df$Cage,
    Regimen = case_when(
      grepl("EB|EBO", Cage) ~ "O-type",
      grepl("CB|CBO", Cage) ~ "B-type",
      TRUE ~ NA_character_
    ),
    Ancestry = case_when(
      grepl("CBO|EBO", Cage) ~ "BO",
      grepl("CB|EB", Cage) ~ "B",
      TRUE ~ NA_character_
    ),
    Sex = summary_df$Sex,
    Replicate = as.factor(gsub("\\D", "", Cage)),
    Treatment = case_when(
      grepl("EBO", Cage) ~ "OBO",
      grepl("EB", Cage) ~ "OB",
      grepl("CBO", Cage) ~ "nBO",
      grepl("CB", Cage) ~ "nB",
      TRUE ~ NA_character_
    )
  )

createPlot <- function(ancestry = FALSE) {
  p <- ggplot(calc_df, aes(x = Sex, y = mean_longevity, fill = Regimen)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
    theme_bw() +
    labs(
      title = "Longevity",
      x = NULL,
      y = "Longevity (days)"
    ) +
    theme(legend.position = "none")
  
  if (ancestry) {
    p <- p + facet_grid(~Ancestry)
  }
  
  return(p)
}

longevity_boxplots_gen12 <- createPlot()
longevity_boxplots_anc_gen12 <- createPlot(ancestry = TRUE)
