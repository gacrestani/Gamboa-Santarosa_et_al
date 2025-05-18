library(readxl)
library(tidyverse)

file_path <- "data/phenotypic_data/raw/Dry and Wet Weight Gen 22.xlsx"

input_df <- read_excel(file_path, na = "NA") %>%
  mutate(`Sex Replicate` = as.character(`Sex Replicate`))

summary_df <- input_df %>%
  group_by(`Population Replicates`, Sex) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(
    wet_weight = `Tube and wet flies (mg)` - `Tube (mg)`,
    Regimen = case_when(
      grepl("EB|EBO", `Population Replicates`) ~ "O-type",
      grepl("CB|CBO", `Population Replicates`) ~ "B-type",
      TRUE ~ NA_character_
    ),
    Ancestry = case_when(
      grepl("EBO|CBO", `Population Replicates`) ~ "BO",
      grepl("EB|CB", `Population Replicates`) ~ "B",
      TRUE ~ NA_character_
    ),
    Treatment = case_when(
      grepl("EB", `Population Replicates`) ~ "EB",
      grepl("EBO", `Population Replicates`) ~ "EBO",
      grepl("CB", `Population Replicates`) ~ "CB",
      grepl("CBO", `Population Replicates`) ~ "CBO",
      TRUE ~ NA_character_
    ),
    Replicate = as.factor(gsub("\\D", "", `Population Replicates`)),
    Sex = gsub("Females", "F", Sex),
    Sex = gsub("Males", "M", Sex),
    wet_weight = as.numeric(wet_weight),
    Group = interaction(Sex, Regimen, sep = "_")
  )

calc_df <- summary_df %>%
  mutate(
    Regimen = factor(Regimen),
    Sex = factor(Sex)
  )

createPlot <- function(ancestry = FALSE) {
  p <- ggplot(calc_df, aes(x = Sex, y = wet_weight, fill = Regimen)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
    labs(
      title = "Wet weight",
      x = NULL,
      y = "Wet weight (mg)"
    ) +
    theme_bw() +
    theme(
      legend.position = "none"
    )
  
  if (ancestry) {
    p <- p + facet_grid(~ Ancestry)
  }
  
  return(p)
}

weight_boxplots <- createPlot()
weight_boxplots_anc <- createPlot(ancestry = TRUE)
