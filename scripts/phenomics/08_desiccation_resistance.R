library(readxl)
library(tidyverse)

file_path <- "data/phenomics_data/Desiccation Resistance Gen 22.xlsx"

read_transposed_xlsx <- function(file) {
  df <- read_excel(file, col_names = FALSE)
  dft <- as.data.frame(t(df[-1]), stringsAsFactors = FALSE)
  names(dft) <- df[[1]]
  dft <- as.data.frame(lapply(dft, type.convert))
  return(dft)
}

input_df <- read_transposed_xlsx(file_path)
colnames(input_df) <- c("Population", "Vial", seq(1, 35, by = 1))

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
    Regime = case_when(
      grepl("EB|EBO", Population) ~ "O-type",
      grepl("CB|CBO", Population) ~ "B-type",
      TRUE ~ NA_character_
    ),
    Sex = case_when(
      grepl("F", Sex) ~ "Female",
      grepl("M", Sex) ~ "Male"
    ),
    Ancestry = case_when(
      grepl("EBO|CBO", Population) ~ "BO",
      grepl("EB|CB", Population) ~ "B",
      TRUE ~ NA_character_
    ),
    Treatment = case_when(
      grepl("EB[1-9]", Population) ~ "OB",
      grepl("EBO", Population) ~ "OBO",
      grepl("CB[1-9]", Population) ~ "nB",
      grepl("CBO", Population) ~ "nBO",
      TRUE ~ NA_character_
    ),
    Replicate = as.factor(gsub("\\D", "", Population)),
    Sex = factor(Sex),
    Regime = factor(Regime)
  )

calc_df$Treatment <- factor(calc_df$Treatment, levels = c("OBO",
                                                          "OB",
                                                          "nBO",
                                                          "nB"))


# PLOTS

desiccation_boxplots <- ggplot(calc_df, aes(x = Regime, y = avg_survival, fill = Regime)) +
  geom_boxplot(width = 0.6) +
  scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
  scale_y_continuous(limits = c(5, 22), breaks = seq(10, 20, by = 5)) +
  scale_x_discrete(labels = c(expression("B"["1-10"]), expression("O"["1-10"]))) +
  labs(
    title = "Desiccation resistance",
    x = NULL,
    y = "Survival (hours)"
  ) +
  theme_bw() +
  theme(#axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.position = "none") +
  facet_wrap(~Sex) +
  stat_compare_means(comparisons = list(c("B-type", "O-type")), label = "p.signif", method = "t.test")

desiccation_boxplots_anc <- ggplot(calc_df, aes(x = Treatment, y = avg_survival, fill = Regime)) +
  geom_boxplot(width = 0.6, aes(group = Treatment)) +
  scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
  scale_y_continuous(limits = c(5, 22), breaks = seq(10, 20, by = 5)) +
  scale_x_discrete(labels = c(expression("OBO"["1-5"]),
                              expression("OB"["1-5"]),
                              expression("nBO"["1-5"]),
                              expression("nB"["1-5"]))) +
  
  labs(
    title = "Desiccation resistance",
    x = NULL,
    y = "Survival (hours)"
  ) +
  theme_bw() +
  theme(#axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    legend.position = "none") +
  facet_wrap(~Sex) +
  stat_compare_means(comparisons = list(
    c(1, 2),
    c(3, 4)),
    label = "p.format",
    method = "t.test",
    label.y = 20)

desiccation_boxplots_anc


# Statistical analysis
lm_fit_desiccation <- lm(avg_survival ~ Regime * Ancestry + Sex, data = calc_df)
summary(lm_fit_desiccation)

confint(lm_fit_desiccation)