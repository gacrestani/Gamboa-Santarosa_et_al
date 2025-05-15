source("scripts/phenotypes/ancestry/17_instant_mortality_curves.R")
source("scripts/phenotypes/ancestry/10B_longevity_boxplots.R")
source("scripts/phenotypes/ancestry/11_development_time.R")
source("scripts/phenotypes/ancestry/12_fecundity.R")
source("scripts/phenotypes/ancestry/13_immune_defense.R")
source("scripts/phenotypes/ancestry/14_weight.R")
source("scripts/phenotypes/ancestry/15_starvation_resistance.R")
source("scripts/phenotypes/ancestry/16_desiccation_resistance.R")

library(gridExtra)

grid.arrange(mortality_curves, longevity_boxplots,
             devtime_boxplots, fecundity_boxplots,
             immune_defense_boxplots, weight_boxplots,
             starvation_boxplots, desiccation_boxplots,
             nrow = 4)

