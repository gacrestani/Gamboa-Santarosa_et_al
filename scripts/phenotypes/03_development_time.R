library(readxl)
library(tidyverse)

path <- "data/phenotypic_data/Development Time Gen 20.xlsx"

input_df <- read_excel(path, na = "-", skip = 1, col_types = "text") %>%
  select(-Notes)

hours <- as.numeric(read_excel(path, col_names = F)[1,])
hours <- hours[!is.na(hours)]

input_df[] <- lapply(input_df, as.character)

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
    )
  )

calc_df$Treatment <- factor(calc_df$Treatment, levels = c("OBO",
                                                          "OB",
                                                          "nBO",
                                                          "nB"))

# PLOTS

devtime_boxplots <- ggplot(calc_df, aes(x = Regime, y = average, fill = Regime)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
  theme_bw() +
  scale_y_continuous(limits = c(210, 280), breaks = seq(210, 270, by = 20)) +
  scale_x_discrete(labels = c(expression("B"["1-10"]), expression("O"["1-10"]))) +
  labs(
    title = "Development time",
    x = NULL,
    y = "Dev. time (hours)"
  ) +
  theme(#axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    legend.position = "none") +
  facet_wrap(~Sex) +
  stat_compare_means(comparisons = list(c("B-type", "O-type")), label = "p.signif", method = "t.test", label.y = 266)

devtime_boxplots

devtime_boxplots_anc <- ggplot(calc_df, aes(x = Treatment, y = average, fill = Regime)) +
  geom_boxplot(width = 0.6, aes(group = Treatment)) +
  scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
  scale_y_continuous(limits = c(210, 280), breaks = seq(210, 270, by = 20)) +
  scale_x_discrete(labels = c(expression("OBO"["1-5"]),
                              expression("OB"["1-5"]),
                              expression("nBO"["1-5"]),
                              expression("nB"["1-5"]))) +
  theme_bw() +
  labs(
    title = "Development time",
    x = NULL,
    y = "Dev. time (hours)"
  ) +
  theme(#axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    legend.position = "none") +
  facet_wrap(~Sex) +
  stat_compare_means(comparisons = list(
    c(1, 2),
    c(3, 4)),
    label = "p.format",
    method = "t.test",
    label.y = 270)


devtime_boxplots_anc

# Statistican analysis
lm_fit_devtime <- lm(average ~ Regime * Ancestry + Sex , data = calc_df)
summary(lm_fit_devtime)

confint(lm_fit_devtime)

