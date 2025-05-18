library(readxl)
library(tidyverse)

file_path <- "data/phenotypic_data/raw/Immunity Gen 22.xlsx"
sheets <- excel_sheets(file_path)
days <- 1:13

getSubgroup <- function(treatment = "Infected") {
  input_df <- read_excel(file_path, sheet = treatment) %>%
    filter(Sex != "UK") %>%
    mutate(Cage = str_remove(Cage, " \\(.*\\)")) %>%
    group_by(Cage, Sex) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop")
  
  calc_matrix <- input_df[, 3:ncol(input_df)] %>% as.matrix()
  discount <- calc_matrix[, 1] + calc_matrix[, 13]
  death_prop <- (rowSums(calc_matrix) - discount) / rowSums(calc_matrix)
  
  result_df <- tibble(
    death_percentage = as.numeric(death_prop),
    Cage = input_df$Cage,
    Sex = input_df$Sex,
    Regimen = case_when(
      grepl("EB|EBO", Cage) ~ "O-type",
      grepl("CB|CBO", Cage) ~ "B-type",
      TRUE ~ NA_character_
    ),
    Ancestry = case_when(
      grepl("EBO|CBO", Cage) ~ "BO",
      grepl("EB|CB", Cage) ~ "B",
      TRUE ~ NA_character_
    ),
    Replicate = as.factor(gsub("\\D", "", Cage)),
    Treatment = treatment
  )
  
  return(result_df)
}


summary_df <- bind_rows(
  getSubgroup("Control"),
  getSubgroup("Infected")
) %>%
  mutate(
    Sex = factor(Sex),
    Regimen = factor(Regimen),
    Treatment = factor(Treatment),
    Group = paste(Regimen, Treatment, sep = "_")
  )

custom_colors <- c(
  "O-type_Control"   = "#377EB8",
  "O-type_Infected"  = "#377EB8",
  "B-type_Control"   = "#E43A3F",
  "B-type_Infected"  = "#E43A3F"
)

createPlot <- function(ancestry = FALSE) {
  p <- ggplot(summary_df, aes(x = Regimen, y = death_percentage, fill = Group)) +
    geom_boxplot(position = "identity", width = 0.6) +
    geom_hline(yintercept = 0.375, linetype = 2) +
    scale_fill_manual(values = custom_colors) +
    theme_bw() +
    labs(
      title = "Immune defense death percentage",
      x = NULL,
      y = "Death percentage"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  if (ancestry) {
    p <- p + facet_grid(. ~ interaction(Sex, Ancestry))
  } else {
    p <- p +
      facet_wrap(~Sex, strip.position = "bottom") +
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank()
      )
    
  }
  
  return(p)
}

immune_defense_boxplots <- createPlot()
immune_defense_boxplots_anc <- createPlot(ancestry = TRUE)
