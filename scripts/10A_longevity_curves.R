# source("scripts/functions.R")
# 
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tidyverse)
# library(gridExtra)

raw_data <- read_excel("data/phenotypic_data/longevity_gen20-shahrestani-version.xlsx")

data$sign <- NULL
data$Notes <- NULL

days <- as.numeric(colnames(data)[3:ncol(data)])
days <- c(0, days[-1] - days[1])

names(data) <- c("Cage", "Sex", days)
data_calcs <- data[3:ncol(data)]

pop_sizes <- rowSums(data_calcs)

longevity_survivalsum <- t(apply(data_calcs, 1, cumsum))
longevity_survivalratio <- 1 - longevity_survivalsum/pop_sizes

plot_df <- as.data.frame(t(longevity_survivalratio))
rownames(plot_df) <- days

df_long <- plot_df %>%
  pivot_longer(cols = -days, names_to = "pop_cages", values_to = "survival_rate")


#longevity_survivalsum <- t(apply(longevity_calcs, 1, cumsum))
#longevity_survivalratio <- 1 - longevity_survivalsum/pop_sizes


plot_df <- as.data.frame(t(data_calcs))
rownames(plot_df) <- days
colnames(plot_df) <- data$Cage
 
plot_df_female <- plot_df[seq(2,120, by = 2)]
plot_df_male <- plot_df[seq(1,120, by = 2)]
 
plot_df_female$days <- days
plot_df_male$days <- days



days_char <- as.character(days)


# 
# GetSurvivalCurves <- function(df, title, legend = FALSE) {
#   # Reshape the data to a long format

#   
#   # df_long$population <- gsub(" .*", "", df_long$pop_cages)
#   # df_long$population <- gsub("EB", "OB", df_long$population)
#   # df_long$population <- gsub("CB", "nB", df_long$population)
#   # 
#   # df_long$treatment <- gsub("[1-9]", "", df_long$population)
#   # df_long$Replicate <- gsub("[A-Z]", "", df_long$population)
#   # 
#   # df_long$Regimen <- df_long$treatment
#   # df_long$Regimen <- gsub("OBO", "O-type", df_long$Regimen)
#   # df_long$Regimen <- gsub("OB", "O-type", df_long$Regimen)
#   # df_long$Regimen <- gsub("nBO", "B-type", df_long$Regimen)
#   # df_long$Regimen <- gsub("nB", "B-type", df_long$Regimen)
#   #  
#   # # Plot the survival curves
#   # p <- ggplot(df_long, aes(x = days, y = survival_rate, group = pop_cages, color = Regimen)) +
#   #   geom_line(size = 1.2) +
#   #   scale_x_continuous(limits = c(0, 105), breaks = c(0,15,30,45,60,75,90,105)) +
#   #   labs(title = title, x = NULL, y = NULL) +
#   #   theme_minimal() +
#   #   scale_color_manual(values = c("O-type" = "black", "B-type" = "grey"))+
#   #   #theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
#   #   theme(legend.position = "none")
#   # 
#   # return(p)
#   return(df_long)
# }
# 
# male <- GetSurvivalCurves(plot_df_male, title = NULL)
# male
# 
# female <- GetSurvivalCurves(plot_df_female, title = NULL)
# female
# 
# grid.arrange(female, male, nrow=1)
# 
# 
# # ---

days_char <- as.character(days)

df_long <- data %>%
  pivot_longer(cols = -c(Cage, Sex), names_to = "pop_cages", values_to = "survival_rate")



longevity_df <- data.frame(mapply(`*`, data_calcs, days))

data$mean_longevity <- rowSums(longevity_df) / rowSums(data_calcs)






# 
# longevity_male <- longevity[longevity$Sex == "M",]
# longevity_female <- longevity[longevity$Sex == "F",]
# 
# 
# 
# 
# 
# mean_longevity <- rowSums(longevity_calcs * days_count) / rowSums(longevity_calcs)
# 
# longevity_boxplot <- longevity_calcs * days_count
# 
# longevity_means <- rowSums(longevity_boxplot)
# longevity_means <- longevity_means/rowSums(longevity_calcs)
# 
# longevity$population <- gsub(" .*", "", longevity$Cage)
# longevity$population <- gsub("EB", "OB", longevity$population)
# longevity$population <- gsub("CB", "nB", longevity$population)
# longevity$treatment <- gsub("[1-9]", "", longevity$population)
# longevity$Regimen <- longevity$treatment
# longevity$Regimen <- gsub("OBO", "O-type", longevity$Regimen)
# longevity$Regimen <- gsub("OB", "O-type", longevity$Regimen)
# longevity$Regimen <- gsub("nBO", "B-type", longevity$Regimen)
# longevity$Regimen <- gsub("nB", "B-type", longevity$Regimen)
# 
# longevity_male <- longevity[longevity$Sex == "M",]
# longevity_female <- longevity[longevity$Sex == "F",]
# 
# 
# p_longevity_m <- ggplot(longevity_male, aes(x = Regimen, y = mean_longevity, fill = Regimen)) +
#   geom_boxplot() +
#   labs(x = "Regimen", y = "Average Longevity", title = "Longevity by Regimen - Female") +
#   scale_y_continuous(limits = c(20,75), breaks = c(20,30,40,50,60,70,80)) +
#   theme_minimal() +
#   scale_fill_manual(values = c("O-type" = "black", "B-type" = "gray"))
# 
# p_longevity_f <- ggplot(longevity_female, aes(x = Regimen, y = mean_longevity, fill = Regimen)) +
#   geom_boxplot() +
#   labs(x = "Regimen", y = "Average Longevity", title = "Longevity by Regimen - Male") +
#   scale_y_continuous(limits = c(20,75), breaks = c(20,30,40,50,60,70,80)) +
#   theme_minimal() +
#   scale_fill_manual(values = c("O-type" = "black", "B-type" = "gray"))
# 
# p_longevity_m
# p_longevity_f
# 
# grid.arrange(p_longevity_f, p_longevity_m, nrow=1)
# 
# 
# 
# 
# o_male <- plot_df_male %>% dplyr::select(matches("^EB|^days"))
# b_male <- plot_df_male %>% dplyr::select(matches("^CB|^days"))
# 
# p_o_male <- GetSurvivalCurves(o_male, title = NULL)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Male
# obo_male <- plot_df_male %>% dplyr::select(matches("^EBO|^days"))
# ob_male <- plot_df_male %>% select(matches("^EB[1-9]|^days"))
# nbo_male <- plot_df_male %>% select(matches("^CBO|^days"))
# nb_male <- plot_df_male %>% select(matches("^CB[1-9]|^days"))
# 
# p_obo_male <- GetSurvivalCurves(obo_male, title = NULL)
# p_ob_male <- GetSurvivalCurves(ob_male, title = NULL)
# p_nbo_male <- GetSurvivalCurves(nbo_male, title = NULL)
# p_nb_male <- GetSurvivalCurves(nb_male, title = NULL)
# 
# # Female
# obo_female <- plot_df_female %>% select(matches("^EBO|^days"))
# ob_female <- plot_df_female %>% select(matches("^EB[1-9]|^days"))
# nbo_female <- plot_df_female %>% select(matches("^CBO|^days"))
# nb_female <- plot_df_female %>% select(matches("^CB[1-9]|^days"))
# 
# p_obo_female <- GetSurvivalCurves(obo_female, title = NULL)
# p_ob_female <- GetSurvivalCurves(ob_female, title = NULL)
# p_nbo_female <- GetSurvivalCurves(nbo_female, title = NULL)
# p_nb_female <- GetSurvivalCurves(nb_female, title = NULL)
# 
# 
# 
# grid.arrange(p_nbo_female + theme(axis.text.y = element_text()),
#              p_nb_female,
#              p_nbo_male,
#              p_nb_male,
#              # Second row
#              p_obo_female + theme(axis.text.x = element_text(), axis.text.y = element_text()),
#              p_ob_female + theme(axis.text.x = element_text()),
#              p_obo_male + theme(axis.text.x = element_text()),
#              p_ob_male + theme(axis.text.x = element_text()),
#              nrow = 2)
# 
# 
# 
# 
