source("scripts/phenotypes/01_instant_mortality_curves.R")
source("scripts/phenotypes/02_longevity_boxplots.R")
source("scripts/phenotypes/03_development_time.R")
source("scripts/phenotypes/04_fecundity.R")
source("scripts/phenotypes/05_immune_defense.R")
source("scripts/phenotypes/06_weight.R")
source("scripts/phenotypes/07_starvation_resistance.R")
source("scripts/phenotypes/08_desiccation_resistance.R")

library(ggpubr)

default <- ggarrange(mortality_curves,
                     longevity_boxplots,
                     devtime_boxplots,
                     fecundity_boxplots,
                     weight_boxplots,
                     starvation_boxplots,
                     desiccation_boxplots,
                     immune_defense_boxplots,
                     labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
                     common.legend = TRUE,
                     legend = "bottom",
                     ncol = 2, nrow = 4)

ancestry <- ggarrange(mortality_curves_anc,
                      longevity_boxplots_anc,
                      devtime_boxplots_anc,
                      fecundity_boxplots_anc,
                      weight_boxplots_anc,
                      starvation_boxplots_anc,
                      desiccation_boxplots_anc,
                      immune_defense_boxplots_anc,
                      labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
                      common.legend = TRUE,
                      legend = "bottom",
                      ncol = 2, nrow = 4)

ggsave("results/figures/figure2.jpeg",
       plot = default,
       width = 8.5, height = 11, units = "in", dpi = 600)

ggsave("results/figures/supplementary_figure2.jpeg",
       plot = ancestry,
       width = 8.5, height = 11, units = "in", dpi = 600)





source("scripts/phenotypes/09_instant_mortality_curves_gen12.R")
source("scripts/phenotypes/10_longevity_boxplots_gen12.R")
source("scripts/phenotypes/11_development_time_gen14.R")


intermediate <- ggarrange(mortality_curves_gen12,
                            ggarrange(longevity_boxplots_gen12,
                                      devtime_boxplots_gen14,
                                      labels = c("B", "C"),
                                      ncol = 2, nrow = 1),
                          labels = "A",
                          common.legend = TRUE,
                          legend = "bottom",
                          nrow = 2)

intermediate_anc <- ggarrange(mortality_curves_anc_gen12,
                          ggarrange(longevity_boxplots_anc_gen12,
                                    devtime_boxplots_anc_gen14,
                                    labels = c("B", "C"),
                                    ncol = 2, nrow = 1),
                          labels = "A",
                          common.legend = TRUE,
                          legend = "bottom",
                          nrow = 2)

ggsave("results/figures/supplementary_figure1.jpeg",
       plot = intermediate,
       width = 8.5, height = 6, units = "in", dpi = 600)

# ggsave("results/figures/intermediate_phenotypes_anc.png",
#        plot = intermediate_anc,
#        width = 8.5, height = 6, units = "in", dpi = 600)


# STATS

summary(cox_model)
summary(lm_fit_longevity)
summary(lm_fit_devtime)
summary(lm_fit_fecundity)
summary(lm_fit_immune_defense)
summary(lm_fit_weight)
summary(lm_fit_starvation)
summary(lm_fit_desiccation)