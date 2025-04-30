# Phenotype analysis

############### Immunity
# Treatment
data <- read.csv("data/phenotypic_data/immunity_gen22_treat.csv")
data$Cage <- sub("[0-9] .*", "", data$Cage) # Group cages by selection treatment
data <- data[data$Sex == "F", ] # Select only females

data$Total <- rowSums(data[,3:15])
data$percentage <- 1 - (data$X13I / data$Total)
data$percentage <- data$percentage * 100

cb_c <- data[data$Cage=="CB","percentage"]
cbo_c <- data[data$Cage=="CBO","percentage"]
eb_c <- data[data$Cage=="EB","percentage"]
ebo_c <- data[data$Cage=="EBO","percentage"]


# Control
data <- read.csv("data/phenotypic_data/immunity_gen22_control.csv")
data$Cage <- sub("[0-9] .*", "", data$Cage)
data <- data[data$Sex == "F", ]

data$Total <- rowSums(data[,3:15])
data$percentage <- 1 - (data$X13I / data$Total)
data$percentage <- data$percentage * 100

cb_t <- data[data$Cage=="CB","percentage"]
cbo_t <- data[data$Cage=="CBO","percentage"]
eb_t <- data[data$Cage=="EB","percentage"]
ebo_t <- data[data$Cage=="EBO","percentage"]



axis_set <- c("OBO" = 1.5,
              "OB" = 3.5,
              "nBO" = 5.5,
              "nB" = 7.5)

# Background rectangles floats
xleft_start <- 0.5
ybottom_start <- -4
xright_start <- 2.5
ytop_start <- 103
boxplot_step <- 2

# Plot
#par(mar=c(5, 4, 4, 8), xpd=TRUE)
boxplot(ebo_t, ebo_c, eb_t, eb_c, cbo_t, cbo_c, cb_t, cb_c,
        col=rep(c("yellow", "green"), 4),
        ylab = "Mean Death (%)",
        xaxt = "n",
        frame.plot = TRUE)
axis(1, at=axis_set, labels=names(axis_set), side = 1)

# EBO rectangle
rect(xleft = xleft_start,
     ybottom = ybottom_start,
     xright = xright_start,
     ytop = ytop_start,
     col = adjustcolor("darkblue", alpha.f = 0.3),
     border = NA)

# EB rectangle
rect(xleft = 2.5,
     ybottom = ybottom_start,
     xright = xright_start + boxplot_step,
     ytop = ytop_start,
     col = adjustcolor("blue", alpha.f = 0.2),
     border = NA)

# CBO
rect(xleft = 2.5 + boxplot_step,
     ybottom = ybottom_start,
     xright = xright_start + 2*boxplot_step,
     ytop = ytop_start,
     col = adjustcolor("darkred", alpha.f = 0.4),
     border = NA)

# CB
rect(xleft = 2.5 + 2*boxplot_step,
     ybottom = ybottom_start,
     xright = xright_start + 3*boxplot_step,
     ytop = ytop_start,
     col = adjustcolor("red", alpha.f = 0.2),
     border = NA)

legend("right",
       inset= c(-0.3,0),
       legend=c("Control", "Infected"),
       col = c("yellow", "green"),
       pch = 15,
       cex = 1,
       box.lty = 0)

############### Longevity

library(readxl)

# Import dataset - this has to be done because of the alpha and beta characters that are not kept in a csv file. We have to use xlsx and the readxl libraries for this
data <- read_excel("data/phenotypic_data/longevity_gen20-shahrestani-version.xlsx")
data$Notes <- NULL
data <- data[data$Sex == "F", ]
data$Sex <- NULL
data$sign <- NULL

days <- read.csv("data/phenotypic_data/days.csv", header = FALSE)
days <- days$V1

getCounts <- function(data, days) {
  counts <- data.frame(Population = character(0), Obs = numeric(0))
  #list <- c()
  
  for (row in 1:nrow(data)) {
    j <- rep(days, data[row,2:ncol(data)])
    
    for (i in 1:length(j)) {
      new_row <- c(Population = data[row,"Cage"], Obs = j[i])
      counts <- rbind(counts, new_row)
    }
  }
  return(counts)
}

couns <- getCounts(data, days)
colnames(couns) <- c("Population", "Obs")

couns$pop <- sub("[0-9] .*", "", couns$Population)
couns$rep <- as.numeric(gsub(".*?([0-9]+).*", "\\1", couns$Population))
couns$cage <- sub(".* ", "", couns$Population)
couns$ancestry_bo <- sub("CBO|EBO", 1, sub("CB|EB", 0, couns$pop))

cb_counts <- couns[couns$pop == "CB",]
cb_rep1 <- cb_counts[cb_counts$rep == 1,]
cb_rep2 <- cb_counts[cb_counts$rep == 2,]
cb_rep3 <- cb_counts[cb_counts$rep == 3,]
cb_rep4 <- cb_counts[cb_counts$rep == 4,]
cb_rep5 <- cb_counts[cb_counts$rep == 5,]

a <- aggregate(cb_rep1$Obs, by = list(cb_rep1$cage), FUN=mean)
mean(a$x)

# long_aov <- aov(Obs ~ as.factor(pop) + as.factor(ancestry_bo), data = couns)
# summary(long_aov)



popi <- "CB"

get_pop_means <- function(popi) {
  pop_counts <- couns[couns$pop == popi, ]
  
  pop_rep1 <- pop_counts[pop_counts$rep == 1,]
  pop_rep2 <- pop_counts[pop_counts$rep == 2,]
  pop_rep3 <- pop_counts[pop_counts$rep == 3,]
  pop_rep4 <- pop_counts[pop_counts$rep == 4,]
  pop_rep5 <- pop_counts[pop_counts$rep == 5,]
  
  rep1 <- aggregate(pop_rep1$Obs, list(pop_rep1$cage), FUN=mean)
  rep2 <- aggregate(pop_rep2$Obs, list(pop_rep2$cage), FUN=mean)
  rep3 <- aggregate(pop_rep3$Obs, list(pop_rep3$cage), FUN=mean)
  rep4 <- aggregate(pop_rep4$Obs, list(pop_rep4$cage), FUN=mean)
  rep5 <- aggregate(pop_rep5$Obs, list(pop_rep5$cage), FUN=mean)
  
  means <- list()
  means <- append(means, mean(rep1$x))
  means <- append(means, mean(rep2$x))
  means <- append(means, mean(rep3$x))
  means <- append(means, mean(rep4$x))
  means <- append(means, mean(rep5$x))
  
  return(means)
}

cb <- unlist(get_pop_means("CB"))
eb <- unlist(get_pop_means("EB"))
ebo <- unlist(get_pop_means("EBO"))
cbo <- unlist(get_pop_means("CBO"))


# ANOVA Analysis
anova_df <- as.data.frame("obs" = c(ebo,eb,cbo,cb),
                          "pop" = c(rep("ebo", 5), rep("eb", 5), rep("cbo", 5), rep("cb", 5))
)

obs <- c(ebo,eb,cbo,cb)
pop <- c(rep("EBO", 5), rep("EB", 5), rep("CBO", 5), rep("CB", 5))
ancestry_bo <- c(rep("BO", 5), rep("B", 5), rep("BO", 5), rep("B", 5))

anova_df <- as.data.frame(cbind(obs, pop, ancestry_bo))

summary(aov(obs ~ pop + ancestry_bo, anova_df))


anova_df <- data.frame(
  obs = c(ebo,eb,cbo,cb),
  pop = factor(c(rep("EBO", 5), rep("EB", 5), rep("CBO", 5), rep("CB", 5))),
  ancestry = factor(c(rep("BO", 5), rep("B", 5), rep("BO", 5), rep("B", 5)))
)

factor(anova_df$ancestry)

model <- aov(obs ~ pop + ancestry, data = anova_df)
summary(model)


# No axis
boxplot(ebo, eb, cbo, cb,
        col=c("darkblue", "blue", "darkred", "red"),
        ylab = "Longevity (days)",
        xaxt = "n")

boxplot(ebo, eb, cbo, cb, names = c("EBO", "EB", "CBO", "CB"), col=c("darkblue", "blue", "darkred", "red"), ylab = "Longevity (days)")


