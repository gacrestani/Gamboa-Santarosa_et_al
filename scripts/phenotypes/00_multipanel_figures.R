source("scripts/phenotypes/instant_mortality_curves.R")
source("scripts/phenotypes/longevity_boxplots.R")
source("scripts/phenotypes/development_time.R")
source("scripts/phenotypes/fecundity.R")
source("scripts/phenotypes/immune_defense.R")
source("scripts/phenotypes/weight.R")
source("scripts/phenotypes/starvation.R")
source("scripts/phenotypes/desiccation.R")

library(ggpubr)

default <- ggarrange(mortality_curves,
                     longevity_boxplots,
                     devtime_boxplots,
                     fecundity_boxplots,
                     immune_defense_boxplots,
                     weight_boxplots,
                     starvation_boxplots,
                     desiccation_boxplots,
                     labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
                     common.legend = TRUE,
                     legend = "bottom",
                     ncol = 2, nrow = 4)

ancestry <- ggarrange(mortality_curves_anc,
                      longevity_boxplots_anc,
                      devtime_boxplots_anc,
                      fecundity_boxplots_anc,
                      immune_defense_boxplots_anc,
                      weight_boxplots_anc,
                      starvation_boxplots_anc,
                      desiccation_boxplots_anc,
                      labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
                      common.legend = TRUE,
                      legend = "bottom",
                      ncol = 2, nrow = 4)

# ggsave("/home/crestang/Repositories/genomic_analysis/results/figures/phenotypic_assays/default_phenotype_panels.svg",
#        plot = default,
#        width = 8.5, height = 11, units = "in", dpi = 600)
# 
# ggsave("/home/crestang/Repositories/genomic_analysis/results/figures/phenotypic_assays/ancestry_phenotype_panels.svg",
#        plot = ancestry,
#        width = 8.5, height = 11, units = "in", dpi = 600)

ggsave("results/default_phenotype_panels.png",
       plot = default,
       width = 8.5, height = 11, units = "in", dpi = 600)

ggsave("results/ancestry_phenotype_panels.png",
       plot = ancestry,
       width = 8.5, height = 11, units = "in", dpi = 600)





source("scripts/phenotypes/instant_mortality_curves_gen12.R")
source("scripts/phenotypes/longevity_boxplots_gen12.R")
source("scripts/phenotypes/development_time_gen14.R")


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

# ggsave("/home/crestang/Repositories/genomic_analysis/results/figures/phenotypic_assays/intermediate_phenotypes.svg",
#        plot = intermediate,
#        width = 8.5, height = 6, units = "in", dpi = 600)
# 
# ggsave("/home/crestang/Repositories/genomic_analysis/results/figures/phenotypic_assays/intermediate_phenotypes_anc.svg",
#        plot = intermediate_anc,
#        width = 8.5, height = 6, units = "in", dpi = 600)

ggsave("results/intermediate_phenotypes.png",
       plot = intermediate,
       width = 8.5, height = 6, units = "in", dpi = 600)

ggsave("results/intermediate_phenotypes_anc.png",
       plot = intermediate_anc,
       width = 8.5, height = 6, units = "in", dpi = 600)
