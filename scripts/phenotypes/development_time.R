library(readxl)
library(tidyverse)

input_df <- read_excel("data/phenotypic_data/raw/Development Time Gen 20.xlsx", skip = 1, col_types = "text") %>%
  select(-Notes)

input_df[] <- lapply(input_df, as.character)

hours <- seq(from = 0, to = 198, by = 6)
colnames(input_df) <- c("Population", "Vials", "Sex", paste0("t", hours))

input_df[input_df == "-"] <- NA

timepoint_cols <- grep("^t", names(input_df), value = TRUE)
input_df[timepoint_cols] <- lapply(input_df[timepoint_cols], as.numeric)

summary_df <- input_df %>%
  group_by(Population, Sex) %>%
  summarise(across(all_of(timepoint_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

calc_matrix <- summary_df[timepoint_cols]
pop_size <- rowSums(calc_matrix, na.rm = TRUE)

weighted_sum <- as.data.frame(mapply(`*`, calc_matrix, hours))
total_time <- rowSums(weighted_sum, na.rm = TRUE)

calc_df <- summary_df %>%
  mutate(
    average = total_time / pop_size,
    Regimen = case_when(
      grepl("EBO|EB", Population) ~ "O-type",
      grepl("CBO|CB", Population) ~ "B-type",
      TRUE ~ NA_character_
    ),
    Ancestry = case_when(
      grepl("EBO|CBO", Population) ~ "BO",
      grepl("EB|CB", Population) ~ "B",
      TRUE ~ NA_character_
    )
  )

createPlot <- function(ancestry = FALSE) {
  p <- ggplot(calc_df, aes(x = Sex, y = average, fill = Regimen)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
    theme_bw() +
    labs(
      title = "Mean development time",
      x = NULL,
      y = "Mean development time (hours)"
    ) +
    theme(legend.position = "none")
  
  if (ancestry) {
    p <- p + facet_grid(~Ancestry)
  }
  
  return(p)
}

devtime_boxplots <- createPlot()
devtime_boxplots_anc <- createPlot(ancestry = TRUE)