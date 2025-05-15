library(readxl)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

file_path <- "data/phenotypic_data/raw/Fecundity Gen 20.xlsx"
sheets <- excel_sheets(file_path)
sheets <- sheets[4:length(sheets)]

# These are your column names in each sheet
vial_cols <- c("Population-Vial #...1", "Population-Vial #...4", "Population-Vial #...7", "Population-Vial #...10")
egg_cols  <- c("# of eggs...2", "# of eggs...5", "# of eggs...8", "# of eggs...11")

needed_cols <- c(vial_cols, egg_cols)

read_and_reshape <- function(sheet_name) {
  df <- read_excel(file_path, sheet = sheet_name, skip = 2)
  
  df <- df %>%
    select(any_of(needed_cols)) %>%
    mutate(across(everything(), as.character))
  
  # Rename columns by position to guarantee correct order
  colnames(df) <- c("vial_1", "vial_2", "vial_3", "vial_4", 
                    "eggs_1", "eggs_2", "eggs_3", "eggs_4")
  
  # Manually reshape: build long dataframe by stacking each pair
  bind_rows(
    df %>% transmute(vial = vial_1, eggs = eggs_1),
    df %>% transmute(vial = vial_2, eggs = eggs_2),
    df %>% transmute(vial = vial_3, eggs = eggs_3),
    df %>% transmute(vial = vial_4, eggs = eggs_4)
  ) %>%
    filter(!(is.na(vial) & is.na(eggs))) %>%
    mutate(sheet = sheet_name)
}

# Apply to all sheets and combine
combined_data <- map_dfr(sheets, read_and_reshape)

combined_data <- combined_data %>%
  separate(vial, into = c("Population", "Vial"), sep = "-")

combined_data <- combined_data %>% mutate(eggs = suppressWarnings(as.numeric(eggs)))

df_avg <- combined_data %>%
  group_by(Population) %>%
  summarise(mean_eggs = mean(eggs, na.rm = TRUE))

df_avg$eggs_average <- df_avg$mean_eggs / 4 # number of females per vial

df_avg <- df_avg %>% mutate(
  Regimen = case_when(
    grepl("EB|EBO", Population) ~ "O-type",
    grepl("CB|CBO", Population) ~ "B-type",
    TRUE ~ NA_character_
  )
)
df_avg$Replicate <- as.factor(gsub("\\D", "", df_avg$Population))
df_avg <- df_avg %>% mutate(
  Treatment = case_when(
    grepl("EB", Population) ~ "EB",
    grepl("EBO", Population) ~ "EBO",
    grepl("CB", Population) ~ "CB",
    grepl("CBO", Population) ~ "CBO",
    TRUE ~ NA_character_
  )
)


fecundity_boxplots <- ggplot(df_avg, aes(x = Regimen, y = eggs_average, fill = Regimen)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  #geom_jitter(width = 0.2, alpha = 0.3, color = "black") +
  scale_fill_manual(values = c("B-type" = "#e43a3f",  # blue
                               "O-type" = "#377EB8")) +# red
  labs(title = "Average number of eggs per female", y = "Average number of eggs per female", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_blank())

fecundity_boxplots
