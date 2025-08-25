source("scripts/phenotypes/instant_mortality_curves.R")
source("scripts/phenotypes/longevity_boxplots.R")
source("scripts/phenotypes/development_time.R")
source("scripts/phenotypes/fecundity.R")
source("scripts/phenotypes/immune_defense.R")
source("scripts/phenotypes/weight.R")
source("scripts/phenotypes/starvation.R")
source("scripts/phenotypes/desiccation.R")

library(gridExtra)

default <- grid.arrange(mortality_curves, longevity_boxplots,
                        devtime_boxplots, fecundity_boxplots,
                        immune_defense_boxplots, weight_boxplots,
                        starvation_boxplots, desiccation_boxplots,
                        nrow = 4)

ancestry <- grid.arrange(mortality_curves_anc, longevity_boxplots_anc,
                         devtime_boxplots_anc, fecundity_boxplots_anc,
                         immune_defense_boxplots_anc, weight_boxplots_anc,
                         starvation_boxplots_anc, desiccation_boxplots_anc,
                         nrow = 4)

ggsave("/var/home/crestang/files/work/research/genomic_analysis/results/figures/phenotypic_assays/default_phenotype_panels.svg",
       plot = default,
       width = 8.5, height = 11, units = "in", dpi = 600)

ggsave("/var/home/crestang/files/work/research/genomic_analysis/results/figures/phenotypic_assays/default_phenotype_panels.png",
       plot = default,
       width = 8.5, height = 11, units = "in", dpi = 600)

ggsave("/var/home/crestang/files/work/research/genomic_analysis/results/figures/phenotypic_assays/ancestry_phenotype_panels.svg",
       plot = ancestry,
       width = 8.5, height = 11, units = "in", dpi = 600)

ggsave("/var/home/crestang/files/work/research/genomic_analysis/results/figures/phenotypic_assays/ancestry_phenotype_panels.png",
       plot = ancestry,
       width = 8.5, height = 11, units = "in", dpi = 600)



source("scripts/phenotypes/instant_mortality_curves_gen12.R")
source("scripts/phenotypes/longevity_boxplots_gen12.R")
source("scripts/phenotypes/development_time_gen14.R")


intermediate <- grid.arrange(mortality_curves_gen12, 
                             longevity_boxplots_gen12, devtime_boxplots_gen14,
                             layout_matrix = rbind(c(1,1),
                                                   c(2,3)))

intermediate_anc <- grid.arrange(mortality_curves_anc_gen12,
                                 longevity_boxplots_anc_gen12, devtime_boxplots_anc_gen14,
                                 layout_matrix = rbind(c(1,1),
                                                       c(2,3)))

ggsave("/var/home/crestang/files/work/research/genomic_analysis/results/figures/phenotypic_assays/intermediate_phenotypes.svg",
       plot = intermediate,
       width = width, height = height, units = "in", dpi = 600)

ggsave("/var/home/crestang/files/work/research/genomic_analysis/results/figures/phenotypic_assays/intermediate_phenotypes_anc.svg",
       plot = intermediate_anc,
       width = width, height = height, units = "in", dpi = 600)
