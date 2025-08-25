library(readxl)
library(tidyverse)

file_path <- "data/phenotypic_data/raw/Starvation Resistance Gen 22.xlsx"

read_transposed_xlsx <- function(file) {
  df <- read_excel(file, col_names = FALSE)
  dft <- as.data.frame(t(df[-1]), stringsAsFactors = FALSE)
  names(dft) <- df[[1]]
  dft <- as.data.frame(lapply(dft, type.convert))
  return(dft)
}

input_df <- read_transposed_xlsx(file_path)
colnames(input_df) <- c("Population", "Vial", seq(from = 3, to = 180, by = 3))

summary_df <- input_df %>%
  separate(Vial, into = c("Sex", "Vial"), sep = 1) %>%
  mutate(Vial = as.character(Vial)) %>%
  select(where(~ !all(is.na(.)))) %>%
  group_by(Population, Sex) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

df_calc <- summary_df[, 3:ncol(summary_df)]
hour_vals <- as.numeric(colnames(summary_df)[3:ncol(summary_df)])

death_counts <- t(apply(df_calc, 1, function(x) {
  x <- as.numeric(x)
  c(x[1], diff(x))
}))

avg_survival <- apply(death_counts, 1, function(deaths) {
  total <- sum(deaths, na.rm = TRUE)
  if (total == 0) return(NA)
  weighted_sum <- sum(deaths * hour_vals, na.rm = TRUE)
  weighted_sum / total
})

calc_df <- summary_df %>%
  mutate(
    avg_survival = avg_survival,
    Regimen = case_when(
      grepl("EB|EBO", Population) ~ "O-type",
      grepl("CB|CBO", Population) ~ "B-type",
      TRUE ~ NA_character_
    ),
    Ancestry = case_when(
      grepl("EBO|CBO", Population) ~ "BO",
      grepl("EB|CB", Population) ~ "B",
      TRUE ~ NA_character_
    ),
    Treatment = case_when(
      grepl("EB", Population) ~ "EB",
      grepl("EBO", Population) ~ "EBO",
      grepl("CB", Population) ~ "CB",
      grepl("CBO", Population) ~ "CBO",
      TRUE ~ NA_character_
    ),
    Replicate = as.factor(gsub("\\D", "", Population)),
    Sex = factor(Sex),
    Regimen = factor(Regimen)
  )

createPlot <- function(ancestry = FALSE) {
  p <- ggplot(calc_df, aes(x = Sex, y = avg_survival, fill = Regimen)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
    labs(
      title = "Starvation resistance",
      x = NULL,
      y = "Survival (hours)"
    ) +
    theme_bw() +
    theme(
      legend.position = "none"
    )
  
  if (ancestry) {
    p <- p + facet_grid(~Ancestry)
  }
  
  return(p)
}

starvation_boxplots <- createPlot()
starvation_boxplots_anc <- createPlot(ancestry = TRUE)
