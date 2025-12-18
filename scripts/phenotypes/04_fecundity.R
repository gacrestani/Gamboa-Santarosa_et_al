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

calc_df <- input_df %>%
  group_by(Population) %>%
  summarise(mean_eggs = mean(eggs, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    eggs_average_per_female = mean_eggs / 4,
    Regime = case_when(
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

fecundity_boxplots <- ggplot(calc_df, aes(x = Regime, y = mean_eggs, fill = Regime)) +
  geom_boxplot(width = 0.6) +
  scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
  scale_y_continuous(limits = c(90, 190), breaks = seq(90, 190, by = 20)) +
  scale_x_discrete(labels = c(expression("B"["1-10"]), expression("O"["1-10"]))) +
  theme_bw() +
  labs(
    title = "Mean population fecundity",
    x = NULL,
    y = "Egg count"
  ) +
  theme(#axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.position = "none") +
  stat_compare_means(comparisons = list(c("B-type", "O-type")), label = "p.signif", method = "t.test", label.y = 170)

fecundity_boxplots

fecundity_boxplots_anc <- ggplot(calc_df, aes(x = Treatment, y = mean_eggs, fill = Regime)) +
  geom_boxplot(width = 0.6, aes(group = Treatment)) +
  scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
  scale_y_continuous(limits = c(90, 190), breaks = seq(90, 190, by = 20)) +
  scale_x_discrete(labels = c(expression("OBO"["1-5"]),
                              expression("OB"["1-5"]),
                              expression("nBO"["1-5"]),
                              expression("nB"["1-5"]))) +
  theme_bw() +
  labs(
    title = "Mean population fecundity",
    x = NULL,
    y = "Egg count"
  ) +
  theme(#axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    legend.position = "none") +
  stat_compare_means(comparisons = list(
    c(1, 2),
    c(3, 4)),
    label = "p.format",
    method = "t.test",
    label.y = 180)

fecundity_boxplots_anc


# Statistical analysis
lm_fit_fecundity <- glm(mean_eggs ~ Regime * Ancestry, family = poisson, data = calc_df)
summary(lm_fit_fecundity)

confint(lm_fit_fecundity)




# 
# # CURVES
# curves_df <- input_df
# curves_df$sheet <- as.numeric(factor(curves_df$sheet, levels = unique(sheets)))
# curves_df <- curves_df %>%
#   mutate(
#     Regime = case_when(
#       grepl("EB|EBO", Population) ~ "O-type",
#       grepl("CB|CBO", Population) ~ "B-type",
#       TRUE ~ NA_character_
#     ),
#     Ancestry = case_when(
#       grepl("EBO|CBO", Population) ~ "BO",
#       grepl("EB|CB", Population) ~ "B",
#       TRUE ~ NA_character_
#     ),
#     Replicate = as.factor(gsub("\\D", "", Population)),
#     Treatment = case_when(
#       grepl("EB[1-9]", Population) ~ "OB",
#       grepl("EBO", Population) ~ "OBO",
#       grepl("CB[1-9]", Population) ~ "nB",
#       grepl("CBO", Population) ~ "nBO",
#       TRUE ~ NA_character_
#     ),
#     eggs_per_female = eggs / 4
#   )
# 
# 
# 
# cuves_plot <- ggplot(curves_df, aes(x = sheet, y = eggs_per_female, group = Population, color = Replicate)) +
#   geom_point(alpha = 0.7) +
#   #scale_color_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
#   #scale_x_continuous(breaks = seq(1, length(sheets), by = 2), labels = seq(1, length(sheets), by = 2)) +
#   labs(
#     title = "Fecundity curves",
#     x = "Day",
#     y = "Egg count"
#   ) +
#   theme_bw() +
#   theme(legend.position = "top") +
#   facet_wrap(. ~ Ancestry + Regime)
# 
