#### code/hgcA_analysis/hgcA_abundance.R ####
# Written by Benjamin D. Peterson


#### Get ready. Get set ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(ggpubr)
library(lubridate)
library(patchwork)
library(readxl)
library(tidyverse)
source("code/HCC_plotting_needs.R")


#### Read in hgcA classification ####
tax.data <- read_xlsx("manuscript/SI/peterson_ISME_SI_tables.xlsx",
                      sheet = "Table S4-hgcA_info") %>%
  rename(seqID = `Sequence ID`,
         manual_classification = `Manual classification`,
         predicted_metabolism = `Predicted metabolism`) %>%
  select(seqID, manual_classification, predicted_metabolism) %>%
  filter(predicted_metabolism != "N/A")


#### Read in data ####
hgcA.data <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv") %>%
  filter(rep) %>%
  left_join(tax.data) %>%
  filter(((RM %in% c(286, 300)) & (year(date) %in% c(2017, 2018))) |
           ((RM %in% c(300, 310)) & (year(date) == 2019))) %>%
  group_by(metagenomeID, date, RM, depth, redoxClassification, predicted_metabolism) %>%
  summarize(coverage = sum(coverage)) %>%
  ungroup() %>%
  mutate(redoxClassification = fct_relevel(redoxClassification,
                                           names(color.vector)[c(1:3, 5)]))


#### Individual group abundance info ####
hgcA.data.individual <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv") %>%
  filter(rep) %>%
  left_join(tax.data)

#### Make color vector ####
# hgcA.manual.taxonomy <- read_xlsx("dataEdited/hgcA_analysis/phylogeny/manual_taxonomy.xlsx",
#                                   sheet = "colors_to_use")
color.vector.metabolism <- c("#61C19A", cb.translator[c("vermillion", "skyblue", "yellow")], "gray50")
names(color.vector.metabolism) <- c("Fermenter", "high redox respiratory", "SRB", "methanogen", "unknown")



#### Relativize data ####
total.data <- hgcA.data %>%
  group_by(metagenomeID) %>%
  summarise(total.coverage = sum(coverage))
hgcA.data <- hgcA.data %>%
  left_join(total.data) %>%
  mutate(rel.coverage = coverage / total.coverage * 100)


#### Summarized data ####
summarized.data <- hgcA.data %>%
  group_by(date, RM, depth, redoxClassification, predicted_metabolism) %>%
  summarise(coverage = sum(coverage)) %>%
  spread(key = predicted_metabolism,
         value = coverage) %>%
  as.data.frame()


#### Plot abundance of functional guilds with hgcA ####
plot.functional.guilds.with.hgcA <- function(redox.state.to.use,
                                             abundance.limits,
                                             legend.location = "none") {
  hgcA.data %>%
    filter(redoxClassification == redox.state.to.use) %>%
    ggplot(aes(x = predicted_metabolism,
               y = coverage)) +
    geom_boxplot(outlier.shape = NA) +
    # geom_line(aes(x = predicted_metabolism,
    #               y = coverage,
    #               group = metagenomeID),
    #           col = "gray60") +
    geom_jitter(aes(col = predicted_metabolism),
                width = 0.1) +
    # geom_point(aes(col = predicted_metabolism)) +
    scale_y_continuous(limits = abundance.limits) +
    scale_color_manual(values = color.vector.metabolism) +
    theme_classic() +
    ylab("Abundance of hgcA+ guild (%)") +
    theme(legend.position = legend.location,
          axis.title.x = element_blank(),
          axis.text.x = element_text(colour = "black",
                                     angle = 45,
                                     vjust = 0.55),
          axis.text.y = element_text(colour = "black"))
}

suboxic.hgcA <- plot.functional.guilds.with.hgcA(redox.state.to.use = "suboxic",
                                                 abundance.limits = c(-0.001, 0.6),
                                                 legend.location = c(0.8, 0.8))
no_nitrate.hgcA <- plot.functional.guilds.with.hgcA(redox.state.to.use = "no_nitrate_no_sulfide",
                                                    abundance.limits = c(-0.001, 0.6))
sulfidic.hgcA <- plot.functional.guilds.with.hgcA(redox.state.to.use = "sulfidic",
                                                  abundance.limits = c(-0.001, 6))

pdf("results/hgcA_analysis/abundance_hgcA_functional_guilds_by_redox.pdf",
    width = 9,
    height = 3.5)
ggarrange(suboxic.hgcA, no_nitrate.hgcA, sulfidic.hgcA,
          nrow = 1, ncol = 3)
dev.off()

#### Plot relative abundance of functional guilds with hgcA ####
# hgcA.data %>%
#   # filter(redoxClassification == redox.state.to.use) %>%
#   ggplot(aes(x = predicted_metabolism,
#              y = rel.coverage)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(aes(col = predicted_metabolism),
#               width = 0.1) +
#   scale_y_continuous(limits = c(-0.01, 100.01)) +
#   scale_color_manual(values = color.vector.metabolism) +
#   facet_wrap(~ redoxClassification,
#              nrow = 1) +
#   theme_classic()
