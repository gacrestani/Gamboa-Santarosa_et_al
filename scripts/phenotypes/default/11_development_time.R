library(readxl)
library(ggplot2)
library(tidyverse)

#library(gridExtra)

data <- read_excel("data/phenotypic_data/raw/Development Time Gen 20.xlsx", skip = 1, col_types = "text")
data$Notes <- NULL
data[] <- sapply(data, as.character)

hours <- seq(from = 0, to = 198, by = 6)
colnames(data) <- c("Population", "Vials", "Sex",
                    paste0("t", hours))

data[data == "-"] <- NA

timepoint_cols <- grep("^t", names(data), value = TRUE)
data[timepoint_cols] <- lapply(data[timepoint_cols], function(x) as.numeric(as.character(x)))

data <- data %>%
  group_by(Population, Sex) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

data_calcs <- data[4:ncol(data)]
data_calcs[] <- sapply(data_calcs, as.numeric)

pop_size <- rowSums(data_calcs, na.rm=T)

data_df <- data.frame(mapply(`*`, data_calcs, hours))
data_sums <- rowSums(data_df, na.rm=T)

data$average <- data_sums/pop_size

data$Regimen <- data$Population
data$Regimen <- gsub("[0-9]+$", "", data$Population)

data$Regimen <- gsub("^EBO$|^EB$", "O-type", data$Regimen)
data$Regimen <- gsub("^CBO$|^CB$", "B-type", data$Regimen)

devtime_boxplots <- ggplot(data, aes(x = Sex, y = average, fill = Regimen)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  #geom_jitter(width = 0.2, alpha = 0.3, color = "black") +
  labs(
    title = "Mean development time",
    x = NULL,
    y = "Mean development time (hours)"
  ) +
  scale_fill_manual(values = c("B-type" = "#E41A1C", "O-type" = "#377EB8")) +
  theme_minimal() +
  theme(legend.position = "none")

devtime_boxplots
