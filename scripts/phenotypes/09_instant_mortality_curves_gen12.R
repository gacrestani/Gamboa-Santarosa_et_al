library(readxl)
library(tidyverse)

input_df <- read_excel("data/phenotypic_data/Longevity Gen 12 (cleaned).xlsx") %>%
  select(-Notes, - `Cage #`)

days <- as.numeric(colnames(input_df)[3:ncol(input_df)])
days <- c(0, days[-1] - days[1])
colnames(input_df) <- c("Cage", "Sex", paste0("d", days))

summary_df <- input_df %>%
  group_by(Cage, Sex) %>%
  summarise(across(where(is.numeric), ~ mean(replace(.x, is.na(.x), 0))), .groups = "drop")

compute_log_mortality <- function(counts) {
  counts <- as.matrix(counts)
  totals <- rowSums(counts, na.rm = TRUE)
  cumulative <- t(apply(counts, 1, cumsum))
  alive_prev <- cbind(totals, totals - cumulative[, -ncol(cumulative)])
  rate <- log(counts / alive_prev)
  rate[!is.finite(rate)] <- NA
  as.data.frame(rate)
}

rate_df <- compute_log_mortality(summary_df[3:ncol(summary_df)]) %>%
  mutate(
    Cage = summary_df$Cage,
    Regime = case_when(
      grepl("EB|EBO", Cage) ~ "O-type",
      grepl("CB|CBO", Cage) ~ "B-type",
      TRUE ~ NA_character_
    ),
    Sex = summary_df$Sex,
    Replicate = as.factor(gsub("\\D", "", Cage)),
    Treatment = case_when(
      grepl("EB", Cage) ~ "EB",
      grepl("EBO", Cage) ~ "EBO",
      grepl("CB", Cage) ~ "CB",
      grepl("CBO", Cage) ~ "CBO",
      TRUE ~ NA_character_
    ),
    Ancestry = case_when(
      grepl("CBO|EBO", Cage) ~ "BO",
      grepl("CB|EB", Cage) ~ "B",
      TRUE ~ NA_character_
    ),
    Sample = rownames(.)
  ) %>%
  mutate(
    Sex = case_when(
      grepl("F", Sex) ~ "Female",
      grepl("M", Sex) ~ "Male")
    )

long_df <- rate_df %>%
  pivot_longer(
    cols = starts_with("d"),
    names_to = "Day",
    values_to = "InstantMortality"
  ) %>%
  mutate(Day = as.numeric(gsub("\\D", "", Day)))


createPlot <- function(ancestry = FALSE) {
  p <- ggplot(long_df, aes(x = Day, y = InstantMortality, group = Sample, color = Regime)) +
    geom_line(alpha = 0.6, linewidth = 0.7) +
    facet_grid(~Sex) +
    scale_color_manual(values = c("O-type" = "#377EB8", "B-type" = "#E43A3F")) +
    theme_bw() +
    labs(
      title = "Age-specific mortality rate - generation 12",
      x = "Day",
      y = "Log mortality rate"
    ) +
    theme(legend.position = "none")
  
  if (ancestry) {
    p <- p + facet_grid(. ~ interaction(Sex, Ancestry))
    
  }
  
  return(p)
  
}

mortality_curves_gen12 <- createPlot()
mortality_curves_anc_gen12 <- createPlot(ancestry = TRUE)

# Statistical analysis
# Cox proportional hazards model
library(survival)
library(survminer)

input_df$N <- rowSums(input_df[3:ncol(input_df)])
input_df$Population <- str_remove(input_df$Cage, " \\(.*\\)")
input_df$Cage <- str_extract(input_df$Cage, "(?<=\\().+?(?=\\))")

input_df <- input_df %>% mutate(
  Regime = case_when(
    grepl("EB|EBO", Population) ~ "O-type",
    grepl("CB|CBO", Population) ~ "B-type",
    TRUE ~ NA_character_
  ),
  Sex = Sex,
  Replicate = as.factor(gsub("\\D", "", Population)),
  Treatment = case_when(
    grepl("EB", Population) ~ "EB",
    grepl("EBO", Population) ~ "EBO",
    grepl("CB", Population) ~ "CB",
    grepl("CBO", Population) ~ "CBO",
    TRUE ~ NA_character_
  ),
  Ancestry = case_when(
    grepl("CBO|EBO", Population) ~ "BO",
    grepl("CB|EB", Population) ~ "B",
    TRUE ~ NA_character_
  )
)

long_input_df <- input_df %>%
  pivot_longer(cols = starts_with("d"), names_to = "Day", values_to = "Dead") %>%
  mutate(Day = as.numeric(str_remove(Day, "d")))

dead_flies <- long_input_df %>%
  filter(Dead > 0) %>%
  uncount(Dead) %>%
  mutate(status = 1) %>% # Dead : status = 1
  select(Cage, Sex, Regime, Ancestry, Replicate, Day, status)


cox_model <- coxph(Surv(Day, status) ~ Regime + Ancestry + Sex, data = dead_flies)
cox.zph(cox_model) # Proportionality is violated

cox_model <- coxph(Surv(Day, status) ~ tt(Regime) + tt(Ancestry) + tt(Sex), data = dead_flies)
summary(cox_model)