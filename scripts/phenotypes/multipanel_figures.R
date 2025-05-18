source("scripts/phenotypes/instant_mortality_curves.R")
source("scripts/phenotypes/longevity_boxplots.R")
source("scripts/phenotypes/development_time.R")
source("scripts/phenotypes/fecundity.R")
source("scripts/phenotypes/immune_defense.R")
source("scripts/phenotypes/weight.R")
source("scripts/phenotypes/starvation.R")
source("scripts/phenotypes/desiccation.R")

library(gridExtra)

grid.arrange(mortality_curves, longevity_boxplots,
             devtime_boxplots, fecundity_boxplots,
             immune_defense_boxplots, weight_boxplots,
             starvation_boxplots, desiccation_boxplots,
             nrow = 4)

grid.arrange(mortality_curves_anc, longevity_boxplots_anc,
             devtime_boxplots_anc, fecundity_boxplots_anc,
             immune_defense_boxplots_anc, weight_boxplots_anc,
             starvation_boxplots_anc, desiccation_boxplots_anc,
             nrow = 4)
