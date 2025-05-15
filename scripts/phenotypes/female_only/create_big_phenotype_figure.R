source("scripts/17_instant_mortality_curves.R")
source("scripts/10B_longevity_boxplots.R")
source("scripts/11_development_time.R")
source("scripts/12_fecundity.R")
source("scripts/13_immune_defense.R")
source("scripts/14_weight.R")
source("scripts/15_starvation_resistance.R")
source("scripts/16_desiccation_resistance.R")

library(gridExtra)

grid.arrange(mortality_curves, longevity_boxplots,
             devtime_boxplots, fecundity_boxplots,
             immune_defense_boxplots, weight_boxplots,
             starvation_boxplots, desiccation_boxplots,
             nrow = 4)

