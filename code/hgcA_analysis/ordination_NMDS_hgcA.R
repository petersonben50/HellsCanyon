#### code/hgcA_analysis/ordination_NMDS_hgcA.R ####
# Written by Benjamin D. Peterson


#### Get ready. Get set ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(lubridate)
library(phyloseq)
library(readxl)
library(tidyverse)
# cb.translator <- readRDS("references/colorblind_friendly_colors.rds")
source("code/HCC_plotting_needs.R")
color.vector <- color.vector[c(2, 3, 5)]
# shape.vector <- c(16, 1, 17, 2, 15, 0)
# names(shape.vector) <- c("2017-RM286", "2017-RM300",
#                          "2018-RM286", "2018-RM300",
#                          "2019-RM300", "2019-RM310")


#### Read in hgcA classification ####
tax.data <- read_xlsx("dataEdited/hgcA_analysis/hgcA_information_edited.xlsx") %>%
  select(seqID, clstr, usedForAbundance, predicted_metabolism)


#### Read in abundance data ####
hgcA.data <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv") %>%
  left_join(tax.data) %>%
  filter(usedForAbundance,
         ((RM %in% c(286, 300)) & (year(date) %in% c(2017, 2018))) |
           ((RM %in% c(300, 310)) & (year(date) == 2019)),
         redoxClassification != "oxic") %>%
  mutate(year = as.character(year(date)))



#### Prepare data for ordination ####
hgcA.data.for.ordination <- hgcA.data %>%
  select(clstr, metagenomeID, coverage) %>%
  spread(key = clstr,
         value = coverage,
         fill = 0)
# Generate phyloseq object
row.names(hgcA.data.for.ordination) <- hgcA.data.for.ordination$metagenomeID
hgcA.data.for.ordination <- hgcA.data.for.ordination %>%
  select(-metagenomeID)
hgcA.phyloseq <- otu_table(hgcA.data.for.ordination,
                           taxa_are_rows = FALSE)


#### Make metagenome metadata vector ####
metadata.vector <- paste(year(hgcA.data$date), "\n",
                         "RM", hgcA.data$RM, ",", hgcA.data$depth, "m",
                         sep = '')
names(metadata.vector) <- hgcA.data$metagenomeID


#### Make metadata a phyloseq class object ####
hgcA.metadata <- hgcA.data %>%
  select(metagenomeID, year, RM, depth, redoxClassification) %>%
  mutate(year.RM = paste(year, "-RM", RM,
                         sep = "")) %>%
  unique()
rownames(hgcA.metadata) <- hgcA.metadata[, "metagenomeID"] %>%
  unlist()
meta.phylo <- sample_data(hgcA.metadata %>% select(-metagenomeID))


#### Combine phyloseq objects ####
all.phyloseq <- merge_phyloseq(hgcA.phyloseq, meta.phylo)


#### Generate ordination ####
hgcA.bray <- ordinate(all.phyloseq,
                      method = "PCoA",
                      distance = "bray")

# bc.ord.hgcA <- 
  plot_ordination(all.phyloseq,
                               hgcA.bray,
                               color = "redoxClassification",
                               shape = "year") +
  scale_color_manual(values = color.vector,
                     labels = renaming.vector,
                     name = "Redox status") +
  scale_shape_manual(values = shape.vector,
                     name = "Year") +
  geom_point(size=5) +
  ylim(c(-0.5, 0.4)) +
  xlim(c(-0.5, 0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        # legend.position = c(0.6, 0.85)
        )

