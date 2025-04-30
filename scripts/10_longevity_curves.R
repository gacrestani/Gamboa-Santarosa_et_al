source("scripts/functions.R")

library("readxl")
library("tidyr")
library("dplyr")
library(gridExtra)

longevity <- readxl::read_xlsx("data/phenotypic_data/longevity_gen20-shahrestani-version.xlsx")

longevity$sign <- NULL
longevity$Notes <- NULL

# Get day_count
days <- as.numeric(colnames(longevity)[3:ncol(longevity)])
days_count <- c(0, days[-1] - days[1])

days
days_count

#
longevity_calcs <- longevity[3:ncol(longevity)]
colnames(longevity_calcs) <- days_count

pop_sizes <- rowSums(longevity_calcs)

longevity_survivalsum <- t(apply(longevity_calcs, 1, cumsum))
longevity_survivalratio <- 1 - longevity_survivalsum/pop_sizes

plot_df <- as.data.frame(t(longevity_survivalratio))
rownames(plot_df) <- days_count

plot_df_female <- plot_df[seq(2,120, by = 2)]
plot_df_male <- plot_df[seq(1,120, by = 2)]

colnames(plot_df_female) <- longevity$Cage[longevity$Sex == "F"]
colnames(plot_df_male) <- longevity$Cage[longevity$Sex == "M"]

plot_df_female$days <- days_count
plot_df_male$days <- days_count


GetSurvivalCurves <- function(df, title, legend = FALSE) {
  # Reshape the data to a long format
  df_long <- df %>%
    pivot_longer(cols = -days, names_to = "pop_cages", values_to = "survival_rate")
  
  df_long$population <- gsub(" .*", "", df_long$pop_cages)
  df_long$population <- gsub("EB", "OB", df_long$population)
  df_long$population <- gsub("CB", "nB", df_long$population)
  
  df_long$treatment <- gsub("[1-9]", "", df_long$population)
  df_long$Replicate <- gsub("[A-Z]", "", df_long$population)
  
  # df_long$Regimen <- df_long$treatment
  # df_long$Regimen <- gsub("OBO", "O-type", df_long$Regimen)
  # df_long$Regimen <- gsub("OB", "O-type", df_long$Regimen)
  # df_long$Regimen <- gsub("nBO", "B-type", df_long$Regimen)
  # df_long$Regimen <- gsub("nB", "B-type", df_long$Regimen)
   
  # Plot the survival curves
  p <- ggplot(df_long, aes(x = days, y = survival_rate, group = pop_cages, color = Replicate)) +
    geom_line(size = 1.2) +
    scale_x_continuous(limits = c(0, 100)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
    theme(legend.position = "none")
  
  return(p)
}

# female <- GetSurvivalCurves(plot_df_female, "Female Survival Curves by Regimen")
# male <- GetSurvivalCurves(plot_df_male, "Male Survival Curves by Regimen")
# 
# obo_male <- plot_df_male %>% select(matches("^CBO|^days"))
# obo_male <- GetSurvivalCurves(obo_male, "")
# 
# obo_male
# 
# ggsave(
#   filename = "results/figures/phenotypes/longevity/female.png",
#   plot = female,
#   width = 2160,
#   height = 1440,
#   bg = "white",
#   units = "px"
# )
# 
# ggsave(
#   filename = "results/figures/phenotypes/longevity/male.png",
#   plot = male,
#   width = 2160,
#   height = 1440,
#   bg = "white",
#   units = "px"
# )

# Male
obo_male <- plot_df_male %>% select(matches("^EBO|^days"))
ob_male <- plot_df_male %>% select(matches("^EB[1-9]|^days"))
nbo_male <- plot_df_male %>% select(matches("^CBO|^days"))
nb_male <- plot_df_male %>% select(matches("^CB[1-9]|^days"))

p_obo_male <- GetSurvivalCurves(obo_male, title = NULL)
p_ob_male <- GetSurvivalCurves(ob_male, title = NULL)
p_nbo_male <- GetSurvivalCurves(nbo_male, title = NULL)
p_nb_male <- GetSurvivalCurves(nb_male, title = NULL)

# Female
obo_female <- plot_df_female %>% select(matches("^EBO|^days"))
ob_female <- plot_df_female %>% select(matches("^EB[1-9]|^days"))
nbo_female <- plot_df_female %>% select(matches("^CBO|^days"))
nb_female <- plot_df_female %>% select(matches("^CB[1-9]|^days"))

p_obo_female <- GetSurvivalCurves(obo_female, title = NULL)
p_ob_female <- GetSurvivalCurves(ob_female, title = NULL)
p_nbo_female <- GetSurvivalCurves(nbo_female, title = NULL)
p_nb_female <- GetSurvivalCurves(nb_female, title = NULL)



grid.arrange(p_nbo_female + theme(axis.text.y = element_text()),
             p_nb_female,
             p_nbo_male,
             p_nb_male,
             # Second row
             p_obo_female + theme(axis.text.x = element_text(), axis.text.y = element_text()),
             p_ob_female + theme(axis.text.x = element_text()),
             p_obo_male + theme(axis.text.x = element_text()),
             p_ob_male + theme(axis.text.x = element_text()),
             nrow = 2)