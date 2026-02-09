library(readxl)
library(tidyverse)

input_df <- read_excel("data/phenotypic_data/Development Time Gen 14.xlsx", skip = 1, col_types = "text", na = "NA")

time <- read_excel("data/phenotypic_data/Development Time Gen 14.xlsx")[1,]
hours <- colnames(time)
hours <- as.numeric(hours[4:length(hours)])

#input_df[] <- lapply(input_df, as.character)

colnames(input_df) <- c("Population", "Vials", "Sex", paste0("t", hours))

timepoint_cols <- grep("^t", names(input_df), value = TRUE)
input_df[timepoint_cols] <- lapply(input_df[timepoint_cols], as.numeric)
input_df$pop_size <- rowSums(input_df[timepoint_cols], na.rm = TRUE)

summary_df <- input_df %>%
  group_by(Population, Sex) %>%
  summarise(across(all_of(c(timepoint_cols, "pop_size")), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
summary_df[summary_df == "NaN"] <- NA


calc_matrix <- summary_df[timepoint_cols]
calc_matrix[calc_matrix == "NaN"] <- NA
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

devtime_boxplots_gen14 <- ggplot(calc_df, aes(x = Regime, y = average, fill = Regime)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
  theme_bw() +
  scale_y_continuous(limits = c(180, 260), breaks = seq(180, 260, by = 20)) +
  scale_x_discrete(labels = c(expression("B"["1-10"]), expression("O"["1-10"]))) +
  labs(
    title = "Development time - generation 14",
    x = NULL,
    y = "Dev. time (hours)"
  ) +
  theme(#axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    legend.position = "none") +
  facet_wrap(~Sex) +
  stat_compare_means(comparisons = list(c("B-type", "O-type")), label = "p.format", method = "t.test", label.y = 240)

devtime_boxplots_anc_gen14 <- ggplot(calc_df, aes(x = Treatment, y = average, fill = Regime)) +
  geom_boxplot(width = 0.6, aes(group = Treatment)) +
  scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
  scale_y_continuous(
    limits = c(180, 260),
    breaks = seq(180, 260, by = 20)) +
  scale_x_discrete(labels = c(expression("OBO"["1-5"]),
                              expression("OB"["1-5"]),
                              expression("nBO"["1-5"]),
                              expression("nB"["1-5"]))) +
  theme_bw() +
  labs(
    title = "Development time - generation 14",
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
    label.y = 240)


# Statistican analysis
lm_fit_devtime <- lm(average ~ Regime * Ancestry + Sex , data = calc_df)
summary(lm_fit_devtime)
