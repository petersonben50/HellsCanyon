####-------------code/binning/prolixibacteraceae_details/abundance_prolix_bins.R-------------####
# Benjamin D. Peterson

# This script will plot the abundance of all the
# Prolixibacteraceae bins from this study.


####-------------Always start with a clean slate-------------####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon/")
library(lubridate)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")
cb.translator <- c(cb.translator, "grey50")
names(cb.translator)[length(cb.translator)] <- "grey"


####-------------Read in depth data-------------####
depth.data <- readRDS("dataEdited/bins/binAnalysis/phylogeny/prolixibacteraceae/depth_of_bins.rds")


####-------------Add in metadata-------------####
MG.metadata <- read_xlsx("metadata/metagenome_metadata.xlsx")
all.data <- depth.data%>%
  left_join(MG.metadata)


####-------------Leaving out HC18ME02_bin_0076 because it's the same as the anvio_hgcA_0130-------------####
all.data <- all.data %>%
  filter(HMS != "HC18ME02_bin_0076")


####-------------Set up color vector-------------####
color.vector <- cb.translator
names(color.vector) <- unique(all.data$binID)


####-------------Plot out data-------------####
all.data %>%
  filter(((RM %in% c(286, 300)) & (year(date) %in% c(2017, 2018))) |
           ((RM %in% c(300, 310)) & (year(date) == 2019)))%>%
  ggplot(aes(y = coverage,
             x = depth,
             group = binID)) +
  geom_point(aes(color = binID)) +
  geom_line(aes(color = binID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~year(date) + RM, nrow = 3) +
  coord_flip(xlim = c(80, 0)) +
  theme_bw()


####-------------Plot out data just from 2018-------------####
all.data %>%
  filter((RM %in% c(286, 300)) & (year(date) == 2018))%>%
  ggplot(aes(y = coverage,
             x = depth,
             group = binID)) +
  geom_point(aes(color = binID)) +
  geom_line(aes(color = binID)) +
  scale_color_manual(values = color.vector) +
  facet_wrap(~year(date) + RM, nrow = 3) +
  coord_flip(xlim = c(80, 0)) +
  theme_bw()
