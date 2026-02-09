library(readxl)
library(tidyverse)
library(ggpubr)

file_path <- "data/phenotypic_data/Immunity Gen 22.xlsx"
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


calc_df <- bind_rows(
  getSubgroup("Control"),
  getSubgroup("Infected")
) %>%
  mutate(
    Sex = case_when(
      grepl("F", Sex) ~ "Female",
      grepl("M", Sex) ~ "Male"
    ),
    Regime = factor(Regime),
    Infection = factor(Treatment),
    Treatment = case_when(
      grepl("EB[1-9]", Cage) ~ "OB",
      grepl("EBO", Cage) ~ "OBO",
      grepl("CB[1-9]", Cage) ~ "nB",
      grepl("CBO", Cage) ~ "nBO",
      TRUE ~ NA_character_
    ),
    Group = paste(Regime, Infection, sep = "_"),
    Group_anc = paste(Group, Ancestry, sep = "_")
  )

calc_df$Treatment <- factor(calc_df$Treatment, levels = c("OBO",
                                                          "OB",
                                                          "nBO",
                                                          "nB"))

calc_df$Group_anc <- factor(calc_df$Group_anc,
                            levels = c("O-type_Control_BO",
                                       "O-type_Control_B",
                                       "O-type_Infected_BO",
                                       "O-type_Infected_B",
                                       "B-type_Control_BO",
                                       "B-type_Control_B",
                                       "B-type_Infected_BO",
                                       "B-type_Infected_B"))



immune_defense_boxplots <- ggplot(calc_df, aes(x = Group, y = death_percentage, fill = Regime)) +
  geom_boxplot(width = 0.6, aes(group = Group), position = "identity") +
  geom_hline(yintercept = 0.5, linetype = 2) +
  scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
  scale_y_continuous(limits = c(0.00, 1.20), breaks = seq(0, 1, by = 0.25)) +
  scale_x_discrete(labels = c("B-type_Control" = expression("B"["1-10"]),
                              "O-type_Control" = expression("O"["1-10"]),
                              "B-type_Infected" = expression("B"["1-10"]),
                              "O-type_Infected" = expression("O"["1-10"]))) +
  theme_bw() +
  labs(
    title = "Immune defense",
    x = NULL,
    y = "Death percentage"
  ) +
  theme(
    legend.position = "none",
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank()
  ) + 
  #annotate("text", x = 1.1, y = 0.425, label = "I.") +
  #annotate("text", x = 1.1, y = 0.325, label = "C.") +
  facet_wrap(~Sex) +
  stat_compare_means(
    comparisons = list(
      c("B-type_Control", "O-type_Control"),
      c("B-type_Infected", "O-type_Infected")),
    label = "p.signif",
    method = "t.test",
    label.y = c(0.2,1))

immune_defense_boxplots

# # Manually calculating p-values so I can manually annotate
# a <- calc_df$death_percentage[
#       calc_df$Regime == "B-type" & 
#       calc_df$Infection == "Infected" &
#       calc_df$Sex == "Female"]
# 
# b <- calc_df$death_percentage[
#       calc_df$Regime == "O-type" & 
#       calc_df$Infection == "Infected" &
#       calc_df$Sex == "Female"]
# 
# inf_f <- t.test(a,b)
# 
# 
# c <- calc_df$death_percentage[
#       calc_df$Regime == "B-type" & 
#       calc_df$Infection == "Control" &
#       calc_df$Sex == "Female"]
# 
# d <- calc_df$death_percentage[
#       calc_df$Regime == "O-type" & 
#       calc_df$Infection == "Control" &
#       calc_df$Sex == "Female"]
# 
# ctrl_f <- t.test(c,d)
# 
# e <- calc_df$death_percentage[
#       calc_df$Regime == "B-type" & 
#       calc_df$Infection == "Infected" &
#       calc_df$Sex == "Male"]
# 
# f <- calc_df$death_percentage[
#       calc_df$Regime == "O-type" & 
#       calc_df$Infection == "Infected" &
#       calc_df$Sex == "Male"]
# 
# inf_m <- t.test(e,f)
# 
# 
# g <- calc_df$death_percentage[
#        calc_df$Regime == "B-type" & 
#         calc_df$Infection == "Control" &
#         calc_df$Sex == "Male"]
# 
# h <- calc_df$death_percentage[
#         calc_df$Regime == "O-type" & 
#         calc_df$Infection == "Control" &
#         calc_df$Sex == "Male"]
# 
# ctrl_m <- t.test(g,h)
# 
# p_text <- data.frame(
#   label = c(paste0("p = ", signif(ctrl_f$p.value, digits = 3)),
#             paste0("p = ", signif(ctrl_m$p.value, digits = 3)),
#             paste0("p = ", signif(inf_f$p.value, digits = 3)),
#             paste0("p = ", signif(inf_m$p.value, digits = 3))),
#   Sex = c("Female", "Male", "Female", "Male"),
#   x = c(1.5, 1.5, 1.5, 1.5),
#   y = c(0.25, 0.25, 1.05, 1.05)
# )
# 
# immune_defense_boxplots <- immune_defense_boxplots +
#   geom_text(
#     data = p_text,
#     mapping = aes(x = x, y = y, label = label)
#   )
# 
# immune_defense_boxplots



immune_defense_boxplots_anc <- ggplot(calc_df, aes(x = Group_anc, y = death_percentage, fill = Regime)) +
  geom_boxplot(width = 0.6, aes(group = Group_anc)) +
  geom_hline(yintercept = 0.375, linetype = 2) +
  scale_fill_manual(values = c("B-type" = "#E43A3F", "O-type" = "#377EB8")) +
  scale_y_continuous(limits = c(0.00, 1.15), breaks = seq(0, 1, by = 0.25)) +
  scale_x_discrete(
    labels = c(expression("OBO"["1-5"]),
               expression("OB"["1-5"]),
               expression("OBO"["1-5"]),
               expression("OB"["1-5"]),
               expression("nBO"["1-5"]),
               expression("nB"["1-5"]),
               expression("nBO"["1-5"]),
               expression("nB"["1-5"]))) +
  theme_bw() +
  labs(
    title = "Immune defense",
    x = NULL,
    y = "Death percentage"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 0.5)
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank()
  ) +
  #facet_wrap(. ~ interaction(Sex, Ancestry))
  #facet_grid(. ~ Sex)
  #annotate("text", x = 1.1, y = 0.425, label = "I.") +
  #annotate("text", x = 1.1, y = 0.325, label = "C.") +
  facet_wrap(~Sex) +
  stat_compare_means(
    comparisons = list(
      c(1, 2),
      c(3, 4),
      c(5, 6),
      c(7, 8)
    ),
    label = "p.format",
    method = "t.test",
    label.y = 1)


immune_defense_boxplots_anc

# Statistical analysis
lm_fit_immune_defense <- glm(cbind(num_dead, num_alive) ~ Regime*Infection + Regime*Ancestry + Sex, family = binomial, data = calc_df)
summary(lm_fit_immune_defense)

confint(lm_fit_immune_defense)
