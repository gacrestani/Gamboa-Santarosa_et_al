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

compute_instantaneous_mortality <- function(dead_mat) {
  # dead_mat is a data frame or matrix with samples as rows, days as columns
  dead_mat <- as.matrix(dead_mat)
  
  # Step 1: total number of flies per sample
  total_flies <- rowSums(dead_mat)
  
  # Step 2: cumulative deaths across days (per row)
  cum_deaths <- t(apply(dead_mat, 1, cumsum))
  
  # Step 3: alive yesterday = total - cumulative deaths up to day-1
  alive_yesterday <- cbind(total_flies, total_flies - cum_deaths[, -ncol(cum_deaths)])
  
  # Step 4: compute ln(dead / alive_yesterday)
  mortality <- log(dead_mat / alive_yesterday)
  
  # Replace infinities and NaNs with NA
  mortality[!is.finite(mortality)] <- NA
  
  
  return(as.data.frame(mortality))
}

calcs <- compute_instantaneous_mortality(calcs)

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




# Step 1: Add a Sample column to mortality_df
calcs$Sample <- rownames(calcs)

# Step 2: Reshape to long format
data_long <- calcs %>%
  pivot_longer(
    cols = starts_with("d"),
    names_to = "Day",
    values_to = "InstantMortality"
  ) %>%
  mutate(
    Day = as.numeric(gsub("\\D", "", Day))  # extract numeric day number
  )

# Step 3: Plot
mortality_curves <- 
  ggplot(data_long, aes(x = Day, y = InstantMortality, group = Sample, color = Regimen)) +
  geom_line(alpha = 0.6, linewidth = 0.7) +
  facet_wrap(~Sex) +
  scale_color_manual(values = c(
      "O-type" = "#377EB8",
      "B-type" = "#E41A1C")) +
  theme_minimal() +
  labs(
    title = "Instantaneous mortality over time",
    x = "Day", y = "Log mortality rate"
  ) +
  theme(legend.position = "none")

mortality_curves
