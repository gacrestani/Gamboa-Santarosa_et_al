library(readxl)

file_path <- "data/phenotypic_data/raw/Desiccation Resistance Gen 22.xlsx"

read.transposed.xlsx <- function(file) {
  df <- read_excel(file, col_names = FALSE)
  dft <- as.data.frame(t(df[-1]), stringsAsFactors = FALSE) 
  names(dft) <- df[,1] 
  dft <- as.data.frame(lapply(dft,type.convert))
  return(dft)            
}

df <- read.transposed.xlsx(file_path)
names(df) <- c("Population", "Vial", seq(from = 1, to = 35, by = 1))

df <- df %>%
  separate(Vial, into = c("Sex", "Vial"), sep = 1)

df$Vial <- as.character(df$Vial)

df <- df[,colSums(is.na(df))<nrow(df)]

df <- df %>%
  group_by(Population, Sex) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

df_calc <- df[,3:ncol(df)]

hour_cols <- colnames(df)[3:ncol(df)]
hour_vals <- as.numeric(hour_cols)

# This will store per-interval deaths (non-cumulative)
death_counts <- t(apply(df_calc, 1, function(x) {
  x <- as.numeric(x)
  c(x[1], diff(x))  # first death count stays, then differences
}))

df$avg_survival <- apply(death_counts, 1, function(deaths) {
  total <- sum(deaths, na.rm = TRUE)
  if (total == 0) return(NA)
  weighted_sum <- sum(deaths * hour_vals, na.rm = TRUE)
  weighted_sum / total
})

df <- df %>% mutate(
  Regimen = case_when(
    grepl("EB|EBO", Population) ~ "O-type",
    grepl("CB|CBO", Population) ~ "B-type",
    TRUE ~ NA_character_
  )
)

df <- df %>% mutate(
  Ancestry = case_when(
    grepl("CBO|EBO", Population) ~ "BO",
    grepl("CB|EB", Population) ~ "B",
    TRUE ~ NA_character_
  )
)

desiccation_boxplots <- ggplot(df, aes(x = Sex, y = avg_survival, fill = Regimen)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  facet_wrap(~ Ancestry) +
  #geom_jitter(width = 0.2, alpha = 0.3, color = "black") +
  labs(
    title = "Desiccation resistance",
    x = NULL,
    y = "Hours survived"
  ) +
  scale_fill_manual(values = c("B-type" = "#E41A1C", "O-type" = "#377EB8")) +
  theme_bw() +
  theme(legend.position = "none")

desiccation_boxplots
