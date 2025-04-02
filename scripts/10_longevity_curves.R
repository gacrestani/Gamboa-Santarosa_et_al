source("scripts/functions.R")

library("readxl")
library("tidyr")
library("dplyr")

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

GetSurvivalCurves <- function(df) {
  # Reshape the data to a long format
  df_long <- df %>%
    pivot_longer(cols = -days, names_to = "pop_cages", values_to = "survival_rate")
  
  df_long$population <- gsub(" .*", "", df_long$pop_cages)
  df_long$population <- gsub("EB", "OB", df_long$population)
  df_long$population <- gsub("CB", "nB", df_long$population)
  
  df_long$treatment <- gsub("[1-9]", "", df_long$population)
  
  df_long$Regimen <- df_long$treatment
  df_long$Regimen <- gsub("OBO", "O-type", df_long$Regimen)
  df_long$Regimen <- gsub("OB", "O-type", df_long$Regimen)
  df_long$Regimen <- gsub("nBO", "B-type", df_long$Regimen)
  df_long$Regimen <- gsub("nB", "B-type", df_long$Regimen)
  
  # Plot the survival curves
  p <- ggplot(df_long, aes(x = days, y = survival_rate, group = pop_cages, color = Regimen)) +
    geom_line(size = 1.2) +
    labs(title = "Female Survival Curves by Regimen", x = "Days", y = "Survival Rate") +
    theme_minimal() + 
    scale_color_manual(values = c("O-type" = "blue",
                                  "B-type" = "red"))
  
  return(p)
}

female <- GetSurvivalCurves(plot_df_female)
male <- GetSurvivalCurves(plot_df_male)

ggsave(
  filename = "results/figures/phenotypes/longevity/female.png",
  plot = female,
  width = 2160,
  height = 1440,
  bg = "white",
  units = "px"
)

ggsave(
  filename = "results/figures/phenotypes/longevity/male.png",
  plot = male,
  width = 2160,
  height = 1440,
  bg = "white",
  units = "px"
)
