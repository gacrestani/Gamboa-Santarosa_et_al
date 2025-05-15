library(readxl)

file_path <- "data/phenotypic_data/raw/Immunity Gen 22.xlsx"
sheets <- excel_sheets(file_path)
days <- 1:13

getSubGroup <- function(treatment = "Infected") {
  df <- read_excel(file_path, sheet = treatment)
  df <- df[which(df$Sex != "UK"),]
  
  df <- df %>%
    mutate(Cage = str_remove(Cage, " \\(.*\\)"))
  
  df <- df %>%
    group_by(Cage, Sex) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop")
  
  df_calc <- df[,3:ncol(df)]
  discount <- df_calc[,1] + df_calc[,13]
  
  result_df <- (rowSums(df_calc) - discount)/rowSums(df_calc)
  names(result_df) <- "death_percentage"
  result_df$Cage <- df$Cage
  result_df$Sex <- df$Sex
  result_df <- result_df %>% mutate(
    Regimen = case_when(
      grepl("EB|EBO", Cage) ~ "O-type",
      grepl("CB|CBO", Cage) ~ "B-type",
      TRUE ~ NA_character_
    )
  )
  result_df$Treatment <- treatment
  
  return(result_df)
}

control <- getSubGroup(treatment = "Control")
infected  <- getSubGroup(treatment = "Infected")

full_df <- rbind(control,infected)

full_df <- full_df %>%
  mutate(
    death_percentage = as.numeric(death_percentage),
    Sex = factor(Sex),
    Regimen = factor(Regimen),
    Treatment = factor(Treatment)
  )

full_df <- full_df %>%
  mutate(Group = paste(Regimen, Treatment, sep = "_"))

custom_colors <- c(
  "O-type_Control"   = "#377EB8",  # light blue
  "O-type_Infected"  = "#377EB8",  # dark blue
  "B-type_Control"   = "#E41A1C",  # light red
  "B-type_Infected"  = "#E41A1C"   # dark red
)


# ggplot(full_df, aes(x = Regimen, y = death_percentage, fill = Group)) +
#   geom_boxplot(position = position_dodge(width = 0.7), width = 0.6) +
#   facet_wrap(~ Sex) +
#   labs(
#     title = "Death % by Regimen, Treatment, and Sex",
#     x = "Regimen",
#     y = "Death Percentage"
#   ) +
#   scale_fill_manual(values = custom_colors, name = "Group") +
#   theme_minimal()

full_df <- full_df %>% mutate(
  Ancestry = case_when(
    grepl("CBO|EBO", Cage) ~ "BO",
    grepl("CB|EB", Cage) ~ "B",
    TRUE ~ NA_character_
  )
)

immune_defense_boxplots <- ggplot(full_df, aes(x = Regimen, y = death_percentage, fill = Group)) +
  geom_boxplot(position = "identity", width = 0.6) +
  #geom_jitter(position = position_jitter(width = 0.2), alpha = 0.3, color = "black") +
  facet_grid(Sex ~ Ancestry) +
  scale_fill_manual(values = custom_colors) +
  #annotate("text", x = 2, y = 0.10, label = "Infected", hjust = 0) +
  labs(
    title = "Immune defense Death percentage",
    x = NULL,
    y = "Death percentage"
  ) +
  geom_hline(yintercept = 0.375, linetype = 2) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank())

immune_defense_boxplots

