source("scripts/functions.R")

library(stringr)
options(scipen=100)

raw_snp_table <- as.data.frame(fread("data/snp_tables/filtered_snps_abcd.txt"))
raw_snp_table[raw_snp_table == "."] <- 0
raw_snp_table[, 6:ncol(raw_snp_table)] <- lapply(raw_snp_table[, 6:ncol(raw_snp_table)], as.numeric)

cov <- GetCov(raw_snp_table)

population <- gsub("N_", "", colnames(cov))

min_cov <- lapply(cov, min)
min_cov <- unlist(min_cov, use.names=FALSE)

mean_cov <- lapply(cov, mean)
mean_cov <- unlist(mean_cov, use.names=FALSE)

reads_df <- read.table("/nfs3/IB/Burke_Lab/Crestani/reads_investigation/read_counts_flylong.txt", header = TRUE, sep = ",")
reads <- as.data.frame(reads_df[grepl("R1", reads_df$Filename),])
reads$population <- str_extract(reads$Filename, "[A-Z]+_rep\\d+_gen\\d+")
reads$Filename <- NULL

library(dplyr)
coverage_df <- reads %>%
  group_by(population) %>% 
  summarise_all(sum)

coverage_df$reads <- coverage_df$ReadCount * 2
coverage_df$min_cov <- min_cov
coverage_df$mean_cov <- mean_cov

coverage_df$expected_cov <- coverage_df$reads * 150 / 120e6

mod <- lm(mean_cov ~ reads, data = coverage_df)
summary(mod)

coverage_df$efficiency <- coverage_df$mean_cov/coverage_df$expected_cov
mean(coverage_df$efficiency)
min(coverage_df$efficiency)
max(coverage_df$efficiency)

#Plot for MSP01
p <- coverage_df |> 
  ggplot(aes(x=reads, y=mean_cov))+ 
  geom_point(size=3, shape=23)+
  scale_y_continuous(name="Mean genome coverage", breaks = seq(0,300,by=50),limits=c(0,300))+
  scale_x_continuous(name="Total reads",labels=scales::comma)+
  theme(
    axis.text.y = element_text(size =12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size=14))+ 
  geom_smooth(method="lm", formula= y~x, col="black", size = 0.5)+
  geom_abline(intercept = 0, slope = (150/180e6), color = "red") +
  annotate("text", x = 200000000, y = 60, label = expression(~italic(p)*"< 2e-16," ~"R^2=0.8844," ~"β=0.0000005"), 
           size = 4, color ="blue")+
  annotate("text", x = 200000000, y = 40, label = "Mean efficiency = 75%", 
           size = 4, color ="blue")+
  ggtitle(label="FlyLong Coverage Report - All Libraries")
p

