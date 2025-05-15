library(readxl)

file_path <- "data/phenotypic_data/raw/Dry and Wet Weight Gen 22.xlsx"
df <- read_excel(file_path, na = "NA")

df$`Sex Replicate` <- as.character(df$`Sex Replicate`)

df <- df %>%
  group_by(`Population Replicates`, Sex) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

df$wet_weight <- df$`Tube and wet flies (mg)` - df$`Tube (mg)`

df <- df %>% mutate(
  Regimen = case_when(
    grepl("EB|EBO", `Population Replicates`) ~ "O-type",
    grepl("CB|CBO", `Population Replicates`) ~ "B-type",
    TRUE ~ NA_character_
  )
)

df <- df %>% mutate(
  Treatment = case_when(
    grepl("EBO", `Population Replicates`) ~ "EBO",
    grepl("EB", `Population Replicates`) ~ "EB",
    grepl("CBO", `Population Replicates`) ~ "CBO",
    grepl("CB", `Population Replicates`) ~ "CB",
    TRUE ~ NA_character_
  )
)
df$Replicate <- as.factor(gsub("\\D", "", df$`Population Replicates`))
df$Sex <- gsub("Females", "F", df$Sex)
df$Sex <- gsub("Males", "M", df$Sex)

df <- df %>%
  mutate(
    wet_weight = as.numeric(wet_weight),
    Regimen = factor(Regimen),
    Sex = factor(Sex),
    Group = interaction(Sex, Regimen, sep = "_")
  )

df <- df %>% mutate(
  Ancestry = case_when(
    grepl("CBO|EBO", Treatment) ~ "BO",
    grepl("CB|EB", Treatment) ~ "B",
    TRUE ~ NA_character_
  )
)

weight_boxplots <- ggplot(df, aes(x = Sex, y = wet_weight, fill = Regimen)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  facet_wrap(~ Ancestry) +
  #geom_jitter(width = 0.2, alpha = 0.3, color = "black") +
  labs(
    title = "Wet weight",
    x = NULL,
    y = "Wet weight (mg)"
  ) +
  scale_fill_manual(values = c("B-type" = "#E41A1C", "O-type" = "#377EB8")) +
  theme_bw() +
  theme(legend.position = "none")

weight_boxplots

