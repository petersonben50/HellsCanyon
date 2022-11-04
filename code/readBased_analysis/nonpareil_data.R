#### code/2019_readBasedAnalyses/nonpareil_data.R ####
# Benjamin D. Peterson

#### Always start with a clean slate ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(gridExtra)
library(Nonpareil)
library(readxl)
library(tidyverse)
cbp.translator <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Read in metadata for metagenomes ####
metadata <- read.csv("metadata/metagenome_metadata.csv") %>%
  mutate(date = mdy(date))


#### Read in sample list  ####
sample_list <- metadata %>%
  mutate(fileLocation = paste("dataEdited/readBasedAnalysis/nonpareil/",
                      metagenomeID,
                      "_NPoutput.npo",
                      sep = "")) %>%
  mutate(locationID = paste(metagenomeID,
                             #" (", year(date), " - ", RM, ":", depth, ")",
                             sep = "")) %>%
  select(fileLocation, locationID)
attach(sample_list)



#### Generate nonpareil object and graph ####
pdf("results/readBasedAnalysis/nonpareil_plots.pdf",
    width = 12,
    height = 6)
nps <- Nonpareil.set(files = fileLocation,
                     labels = locationID)
dev.off()
detach(sample_list)



#### Plot diversity ####
np.diversity <- summary(nps)[,"diversity"]
np.diversity.df <- metadata %>%
  left_join(data.frame(metagenomeID = names(np.diversity),
                       diversity = np.diversity)) %>%
  arrange(depth) %>% arrange(RM) %>% arrange(year(date))

# pdf("results/readBasedAnalysis/nonpareil_diversity.pdf",
#     width = 6,
#     height = 4)
np.diversity.df %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  ggplot(aes(x = siteID,
             y = diversity)) +
  geom_point() +
  theme_classic() +
  ylim(c(20, 25)) +
  ylab("Nonpareil diversity index")
# dev.off()

# Run a one-way ANOVA on diversity
# one.way.anova <- aov(diversity ~ siteID,
#                      data = np.diversity.df)
# sink("results/2019_readBasedAnalysis/nonpareil_diversity_anova.txt")
# summary(one.way.anova)
# sink()
# Not significantly different.



#### Make chart of coverage ####
MGcoverage.vector <- summary(nps)[,"C"]*100
MGcoverage.df <- metadata %>%
  left_join(data.frame(metagenomeID = names(MGcoverage.vector),
                       percent_coverage = MGcoverage.vector)) %>%
  arrange(depth) %>% arrange(RM) %>% arrange(year(date)) %>%
  mutate(percent_coverage = round(percent_coverage, 1))



pdf("results/readBasedAnalysis/nonpareil_coverage.pdf",
    width = 6,
    height = 10)
grid.table(MGcoverage.df,
           rows = NULL,
           cols = colnames(MGcoverage.df))
dev.off()


#### Also save to csv ####
write.csv(MGcoverage.df %>%
            left_join(np.diversity.df %>% 
                        select(metagenomeID, diversity)),
          "results/readBasedAnalysis/nonpareil_coverage.csv",
          row.names = FALSE)
