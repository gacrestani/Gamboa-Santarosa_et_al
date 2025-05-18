library(readxl)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

file_path <- "data/phenotypic_data/raw/Fecundity Gen 20.xlsx"
sheets <- excel_sheets(file_path)[4:length(excel_sheets(file_path))]

vial_cols <- c("Population-Vial #...1", "Population-Vial #...4", "Population-Vial #...7", "Population-Vial #...10")
egg_cols  <- c("# of eggs...2", "# of eggs...5", "# of eggs...8", "# of eggs...11")
needed_cols <- c(vial_cols, egg_cols)

read_and_reshape <- function(sheet_name) {
  df <- read_excel(file_path, sheet = sheet_name, skip = 2) %>%
    select(any_of(needed_cols)) %>%
    mutate(across(everything(), as.character))
  
  colnames(df) <- c("vial_1", "vial_2", "vial_3", "vial_4", 
                    "eggs_1", "eggs_2", "eggs_3", "eggs_4")
  
  bind_rows(
    df %>% transmute(vial = vial_1, eggs = eggs_1),
    df %>% transmute(vial = vial_2, eggs = eggs_2),
    df %>% transmute(vial = vial_3, eggs = eggs_3),
    df %>% transmute(vial = vial_4, eggs = eggs_4)
  ) %>%
    filter(!(is.na(vial) & is.na(eggs))) %>%
    mutate(sheet = sheet_name)
}

input_df <- map_dfr(sheets, read_and_reshape) %>%
  separate(vial, into = c("Population", "Vial"), sep = "-") %>%
  mutate(eggs = suppressWarnings(as.numeric(eggs)))

summary_df <- input_df %>%
  group_by(Population) %>%
  summarise(mean_eggs = mean(eggs, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    eggs_average = mean_eggs / 4,
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
    Replicate = as.factor(gsub("\\D", "", Population)),
    Treatment = case_when(
      grepl("EB", Population) ~ "EB",
      grepl("EBO", Population) ~ "EBO",
      grepl("CB", Population) ~ "CB",
      grepl("CBO", Population) ~ "CBO",
      TRUE ~ NA_character_
    )
  )

calc_df <- summary_df

createPlot <- function(ancestry = FALSE) {
  p <- ggplot(calc_df, aes(x = Regimen, y = eggs_average, fill = Regimen)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
    labs(
      title = "Average number of eggs per female",
      x = NULL,
      y = "Average number of eggs per female"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  if (ancestry) {
    p <- p + facet_grid(~Ancestry)
  }
  
  return(p)
}

fecundity_boxplots <- createPlot()
fecundity_boxplots_anc <- createPlot(ancestry = TRUE)
