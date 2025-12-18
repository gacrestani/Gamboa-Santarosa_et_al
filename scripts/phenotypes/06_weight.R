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
    Regime = case_when(
      grepl("EB|EBO", `Population Replicates`) ~ "O-type",
      grepl("CB|CBO", `Population Replicates`) ~ "B-type",
      TRUE ~ NA_character_
    ),
    Sex = case_when(
      grepl("F", Sex) ~ "Female",
      grepl("M", Sex) ~ "Male"
    ),
    Ancestry = case_when(
      grepl("EBO|CBO", `Population Replicates`) ~ "BO",
      grepl("EB|CB", `Population Replicates`) ~ "B",
      TRUE ~ NA_character_
    ),
    Treatment = case_when(
      grepl("EB[1-9]", `Population Replicates`) ~ "OB",
      grepl("EBO", `Population Replicates`) ~ "OBO",
      grepl("CB[1-9]", `Population Replicates`) ~ "nB",
      grepl("CBO", `Population Replicates`) ~ "nBO",
      TRUE ~ NA_character_
    ),
    Replicate = as.factor(gsub("\\D", "", `Population Replicates`)),
    Sex = gsub("Females", "F", Sex),
    Sex = gsub("Males", "M", Sex),
    wet_weight = as.numeric(wet_weight),
    Group = interaction(Sex, Regime, sep = "_")
  )

calc_df <- summary_df %>%
  mutate(
    Regime = factor(Regime),
    Sex = factor(Sex)
  )

calc_df$Treatment <- factor(calc_df$Treatment, levels = c("OBO",
                                                          "OB",
                                                          "nBO",
                                                          "nB"))



# PLOTS
weight_boxplots <- ggplot(calc_df, aes(x = Regime, y = wet_weight, fill = Regime)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
  scale_y_continuous(limits = c(5, 15), breaks = seq(5, 12.5, by = 2.5)) +
  scale_x_discrete(labels = c(expression("B"["1-10"]), expression("O"["1-10"]))) +
  labs(
    title = "Mean population weight",
    x = NULL,
    y = "Weight (mg)"
  ) +
  theme_bw() +
  theme(#axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.position = "none") +
  facet_wrap(~Sex) +
  stat_compare_means(comparisons = list(c("B-type", "O-type")), label = "p.signif", method = "t.test")


weight_boxplots

weight_boxplots_anc <- ggplot(calc_df, aes(x = Treatment, y = wet_weight, fill = Regime)) +
  geom_boxplot(width = 0.6, aes(group = Treatment)) +
  scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
  scale_y_continuous(limits = c(5, 14), breaks = seq(6, 12, by = 2)) +
  scale_x_discrete(labels = c(expression("OBO"["1-5"]),
                              expression("OB"["1-5"]),
                              expression("nBO"["1-5"]),
                              expression("nB"["1-5"]))) +
  labs(
    title = "Mean population weight",
    x = NULL,
    y = "Weight (mg)"
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
    label.y = 13)


weight_boxplots_anc

# Anova test
lm_fit_weight <- lm(wet_weight ~ Regime * Ancestry + Sex, data = calc_df)
summary(lm_fit_weight)

confint(lm_fit_weight)
