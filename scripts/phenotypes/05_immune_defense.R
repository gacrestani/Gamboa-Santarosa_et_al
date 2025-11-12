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
  
  names(input_df)[length(names(input_df))] <- "num_alive"
  
  calc_matrix <- input_df[, 4:(ncol(input_df)-1)] %>% as.matrix() # First day data censored
  death_prop <- rowSums(calc_matrix) / (rowSums(calc_matrix) + input_df$num_alive)
  
  result_df <- tibble(
    death_percentage = as.numeric(death_prop),
    num_dead = rowSums(calc_matrix),
    num_alive = input_df$num_alive,
    Cage = input_df$Cage,
    Sex = input_df$Sex,
    Regime = case_when(
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
    Regime = factor(Regime),
    Infection = factor(Treatment),
    Group = paste(Regime, Treatment, sep = "_")
  )

custom_colors <- c(
  "O-type_Control"   = "#377EB8",
  "O-type_Infected"  = "#377EB8",
  "B-type_Control"   = "#E43A3F",
  "B-type_Infected"  = "#E43A3F"
)

createPlot <- function(ancestry = FALSE) {
  p <- ggplot(summary_df, aes(x = Regime, y = death_percentage, fill = Group)) +
    geom_boxplot(position = "identity", width = 0.6) + # if you add geom_jitter(), you will have 20 data points per boxplot. That's because of position="identity", and it makes it so that the boxplots better occupy the plot space
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
immune_defense_boxplots +
  annotate("text", x = 0.65, y = 0.425, label = "Infected") +
  annotate("text", x = 0.65, y = 0.325, label = "Control")


immune_defense_boxplots_anc <- createPlot(ancestry = TRUE)

# Statistical analysis
lm_fit_immune_defense <- glm(cbind(num_dead, num_alive) ~ Regime*Infection + Ancestry + Sex, family = binomial, data = summary_df)
summary(lm_fit_immune_defense)

car::Anova(lm_fit_immune_defense)

library(emmeans)
interaction_analysis <- emmeans(lm_fit_immune_defense, ~ Regime | Infection)
interaction_analysis <- emmeans(lm_fit_immune_defense, ~ Regime * Infection, type = "response")
pairs(interaction_analysis, adjust = "tukey")
