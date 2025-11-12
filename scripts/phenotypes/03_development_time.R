library(readxl)
library(tidyverse)

input_df <- read_excel("data/phenotypic_data/raw/Development Time Gen 20.xlsx", skip = 1, col_types = "text") %>%
  select(-Notes)

input_df[] <- lapply(input_df, as.character)

hours <- seq(from = 6 + 189.6, to = 204 + 189.6, by = 6)
colnames(input_df) <- c("Population", "Vials", "Sex", paste0("t", hours))

input_df[input_df == "-"] <- NA

timepoint_cols <- grep("^t", names(input_df), value = TRUE)
input_df[timepoint_cols] <- lapply(input_df[timepoint_cols], as.numeric)
input_df$pop_size <- rowSums(input_df[timepoint_cols], na.rm = TRUE)

summary_df <- input_df %>%
  group_by(Population, Sex) %>%
  summarise(across(all_of(c(timepoint_cols, "pop_size")), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

calc_matrix <- summary_df[timepoint_cols]
pop_size <- summary_df$pop_size

weighted_sum <- as.data.frame(mapply(`*`, calc_matrix, hours))
total_time <- rowSums(weighted_sum, na.rm = TRUE)

calc_df <- summary_df %>%
  mutate(
    average = total_time / pop_size,
    Regime = case_when(
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
  p <- ggplot(calc_df, aes(x = Sex, y = average, fill = Regime)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
    theme_bw() +
    scale_y_continuous(limits = c(210, 270), breaks = seq(210, 270, by = 20)) +
    labs(
      #title = "Mean development time",
      x = NULL,
      y = "Mean development time (hours)"
    ) +
    theme(legend.position = "none")
  
  if (ancestry) {
    p <- p + facet_grid(~Ancestry)
  }
  
  p <- p + stat_compare_means(aes(group = Regime), label = "p.signif", method = "t.test")
  
  return(p)
}

devtime_boxplots <- createPlot()
devtime_boxplots_anc <- createPlot(ancestry = TRUE)

# Statistican analysis
lm_fit_devtime <- lm(average ~ Regime + Ancestry + Sex, data = calc_df)
summary(lm_fit_devtime)
