#### code/hgcA_analysis/hgcA_abundance.R ####
# Written by Benjamin D. Peterson


#### Get ready. Get set ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(lubridate)
library(patchwork)
library(readxl)
library(tidyverse)
source("code/plotting_functions.R")
cb.translator.short <- readRDS("references/colorblind_friendly_colors.rds")
cb.translator <- c(cb.translator.short, "gray50", "gray20", "gray80")
names(cb.translator) <- c(names(cb.translator.short), "gray50", "gray20", "gray80")


#### Read in hgcA classification ####
tax.data <- read_xlsx("dataEdited/hgcA_analysis/hgcA_information_edited.xlsx") %>%
  select(seqID, manual_classification, predicted_metabolism)


#### Make color vector ####
hgcA.manual.taxonomy <- read_xlsx("dataEdited/hgcA_analysis/phylogeny/manual_taxonomy.xlsx",
                                  sheet = "colors_to_use")
color.vector <- cb.translator[hgcA.manual.taxonomy$colorsToUse]
names(color.vector) <- hgcA.manual.taxonomy$seqID



#### Read in data ####
hgcA.data <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv") %>%
  filter(rep) %>%
  left_join(tax.data)



#### Variables to set ####
abundance.cutoff = 0.05
suboxic.fig.width = 4
suboxic.fig.height = 4
noN.noS.width = 2
noN.noS.height = 4
sulf.fig.width = 2
sulf.fig.height = 4


#### Focus on high-hgcA content samples ####
high.hgcA.list <- hgcA.data %>%
  group_by(metagenomeID) %>%
  summarise(coverage = sum(coverage)) %>%
  filter(coverage >= abundance.cutoff)

# #### Set it up so that metabolic potentials are in same order, even if they aren't present in a given sample ####
# hgcA.data <- hgcA.data %>%
#          predicted_metabolism = fct_relevel(predicted_metabolism,
#                                             ))


#### Plot total hgcA content under nitrate-reducing conditions ####
nitrate.total.hgcA <- hgcA.data %>%
  filter(redoxClassification == "suboxic",
         metagenomeID %in% high.hgcA.list$metagenomeID)%>%
  mutate(date_location = paste("RM", RM, "\n", year(date), ", ", depth, "m",
                               sep = "")) %>%
  group_by(date_location) %>%
  summarize(coverage = sum(coverage)) %>%
  ggplot(aes(x = date_location,
             y = coverage)) +
  geom_bar(stat = "identity",
           fill = "white",
           col = "gray25") +
  theme_classic() +
  ylim(0, 1.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf("results/hgcA_analysis/barplots_figure/suboxic_total_hgcA.pdf",
    height = suboxic.fig.height,
    width = suboxic.fig.width)
nitrate.total.hgcA
dev.off()


#### Plot guilds of hgcA content under nitrate-reducing conditions ####
nitrate.types.hgcA <- hgcA.data %>%
  filter(redoxClassification == "suboxic",
         metagenomeID %in% high.hgcA.list$metagenomeID) %>%
  mutate(date_location = paste("RM", RM, "\n", year(date), ", ", depth, "m",
                               sep = "")) %>%
  group_by(date_location, predicted_metabolism, manual_classification) %>%
  summarize(coverage = sum(coverage)) %>%
  # Set it up so that metabolic potentials are in same order, even if they aren't present in that sample.
  spread(key = date_location,
         value = coverage,
         fill = 0) %>%
  gather(key = date_location,
         value = coverage,
         -c(predicted_metabolism, manual_classification)) %>%
  mutate(predicted_metabolism = factor(predicted_metabolism,
                                       levels = c("Fermenter", "high redox respiratory", "SRB", "methanogen", "unknown"))) %>%
  ggplot(aes(x = date_location,
             y = coverage,
             group = predicted_metabolism,
             fill = manual_classification)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  scale_fill_manual(values = color.vector) +
  theme_classic() +
  ylim(0, 1.25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
pdf("results/hgcA_analysis/barplots_figure/suboxic_types_hgcA.pdf",
    height = suboxic.fig.height,
    width = suboxic.fig.width)
nitrate.types.hgcA
dev.off()




#### Plot total of hgcA content under no nitrate, no sulfide conditions ####
no.nitrate.no.sulfide.total.hgcA <- hgcA.data %>%
  filter(redoxClassification == "no_nitrate_no_sulfide",
         metagenomeID %in% high.hgcA.list$metagenomeID) %>%
  mutate(date_location = paste("RM", RM, "\n", year(date), ", ", depth, "m",
                               sep = "")) %>%
  group_by(date_location) %>%
  summarize(coverage = sum(coverage)) %>%
  ggplot(aes(x = date_location,
             y = coverage)) +
  geom_bar(stat = "identity",
           fill = "white",
           col = "gray25") +
  theme_classic() +
  ylim(0, 1.2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf("results/hgcA_analysis/barplots_figure/no_nitrate_no_sulfide_total_hgcA.pdf",
    height = noN.noS.height,
    width = noN.noS.width)
no.nitrate.no.sulfide.total.hgcA
dev.off()


#### Plot guilds of hgcA content under no nitrate, no sulfide conditions ####
no.nitrate.no.sulfide.types.hgcA <- hgcA.data %>%
  filter(redoxClassification == "no_nitrate_no_sulfide",
         metagenomeID %in% high.hgcA.list$metagenomeID) %>%
  mutate(date_location = paste("RM", RM, "\n", year(date), ", ", depth, "m",
                               sep = "")) %>%
  group_by(date_location, predicted_metabolism, manual_classification) %>%
  summarize(coverage = sum(coverage)) %>%
  # Set it up so that metabolic potentials are in same order, even if they aren't present in that sample.
  spread(key = date_location,
         value = coverage,
         fill = 0) %>%
  gather(key = date_location,
         value = coverage,
         -c(predicted_metabolism, manual_classification)) %>%
  mutate(predicted_metabolism = factor(predicted_metabolism,
                                       levels = c("Fermenter", "high redox respiratory", "SRB", "methanogen", "unknown"))) %>%
  ggplot(aes(x = date_location,
             y = coverage,
             group = predicted_metabolism,
             fill = manual_classification)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  scale_fill_manual(values = color.vector) +
  theme_classic() +
  ylim(0, 1.2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
no.nitrate.no.sulfide.types.hgcA
pdf("results/hgcA_analysis/barplots_figure/no_nitrate_no_sulfide_types_hgcA.pdf",
    height = noN.noS.height,
    width = noN.noS.width)
no.nitrate.no.sulfide.types.hgcA
dev.off()





#### Plot total of hgcA content under sulfidic conditions ####
sulfidic.total.hgcA <- hgcA.data %>%
  filter(redoxClassification == "sulfidic",
         metagenomeID %in% high.hgcA.list$metagenomeID) %>%
  mutate(date_location = paste("RM", RM, "\n", year(date), ", ", depth, "m",
                               sep = "")) %>%
  group_by(date_location) %>%
  summarize(coverage = sum(coverage)) %>%
  ggplot(aes(x = date_location,
             y = coverage)) +
  geom_bar(stat = "identity",
           fill = "white",
           col = "gray25") +
  theme_classic() +
  ylim(0, 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf("results/hgcA_analysis/barplots_figure/sulfidic_total_hgcA.pdf",
    height = sulf.fig.height,
    width = sulf.fig.width)
sulfidic.total.hgcA
dev.off()


#### Plot guilds of hgcA content under sulfidic conditions ####
sulfidic.types.hgcA <- hgcA.data %>%
  filter(redoxClassification == "sulfidic",
         metagenomeID %in% high.hgcA.list$metagenomeID) %>%
  mutate(date_location = paste("RM", RM, "\n", year(date), ", ", depth, "m",
                               sep = "")) %>%
  group_by(date_location, predicted_metabolism, manual_classification) %>%
  summarize(coverage = sum(coverage)) %>%
  # Set it up so that metabolic potentials are in same order, even if they aren't present in that sample.
  spread(key = date_location,
         value = coverage,
         fill = 0) %>%
  gather(key = date_location,
         value = coverage,
         -c(predicted_metabolism, manual_classification)) %>%
  mutate(predicted_metabolism = factor(predicted_metabolism,
                                       levels = c("Fermenter", "high redox respiratory", "SRB", "methanogen", "unknown"))) %>%
  ggplot(aes(x = date_location,
             y = coverage,
             group = predicted_metabolism,
             fill = manual_classification)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  scale_fill_manual(values = color.vector) +
  theme_classic() +
  ylim(0, 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
sulfidic.types.hgcA
pdf("results/hgcA_analysis/barplots_figure/sulfidic_types_hgcA.pdf",
    height = sulf.fig.height,
    width = sulf.fig.width)
sulfidic.types.hgcA
dev.off()
