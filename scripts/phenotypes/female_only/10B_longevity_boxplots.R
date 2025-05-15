library(readxl)
library(tidyverse)

raw_data <- read_excel("data/phenotypic_data/longevity_gen20-shahrestani-version.xlsx")

raw_data$sign <- NULL
raw_data$Notes <- NULL

days <- as.numeric(colnames(raw_data)[3:ncol(raw_data)])
days <- c(0, days[-1] - days[1])

names(raw_data) <- c("Cage", "Sex", paste0("d",days))

data <- raw_data %>%
  mutate(Cage = str_remove(Cage, " \\(.*\\)"))

data <- data %>%
  group_by(Cage, Sex) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

calcs <- data[3:ncol(data)]

total_flies <- rowSums(calcs)


longevity <- sweep(calcs, MARGIN=2, days, `*`)
total_longevity <- rowSums(longevity)

calcs$mean_longevity <- total_longevity / total_flies
calcs$Cage <- data$Cage
calcs <- calcs %>% mutate(
  Regimen = case_when(
    grepl("EB|EBO", Cage) ~ "O-type",
    grepl("CB|CBO", Cage) ~ "B-type",
    TRUE ~ NA_character_
  )
)
calcs$Sex <- data$Sex
calcs$Replicate <- as.factor(gsub("\\D", "", calcs$Cage))
calcs <- calcs %>% mutate(
  Treatment = case_when(
    grepl("EB", Cage) ~ "EB",
    grepl("EBO", Cage) ~ "EBO",
    grepl("CB", Cage) ~ "CB",
    grepl("CBO", Cage) ~ "CBO",
    TRUE ~ NA_character_
  )
)

calcs_female <- calcs[which(calcs$Sex == "F"),]

longevity_boxplots <- ggplot(calcs_female, aes(y = mean_longevity, fill = Regimen)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  #geom_jitter(width = 0.2, alpha = 0.3, color = "black") +
  scale_y_continuous(limits = c(30,70)) +
  labs(
    title = "Mean longevity",
    x = NULL,
    y = "Mean longevity (days)"
  ) +
  scale_fill_manual(values = c("B-type" = "#e43a3f", "O-type" = "#377EB8")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_blank())

longevity_boxplots
