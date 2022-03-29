#### code/binning/depth_plots.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
cb.translator <- c(cb.translator, "grey50")
names(cb.translator)[length(cb.translator)] <- "grey"


#### Read in metadata ####
MG.metadata <- read_xlsx("metadata/metagenome_metadata.xlsx")


#### Read in metabolism data ####
metabolic.data <- read_xlsx("dataEdited/binning/metabolism/metabolic_summary.xlsx",
                            sheet = "batch_HMMs") %>%
  select(HMS, binID, metabolic_assignment)
HMS.colors <- read_xlsx("dataEdited/binning/metabolism/metabolic_summary.xlsx",
                        sheet = "HMS_colors") %>%
  select(HMS, colorToUse)
HMS.colors.vector <- cb.translator[HMS.colors$colorToUse]
names(HMS.colors.vector) <- HMS.colors$HMS
rm(HMS.colors)


#### Read in depth data ####
depth.data <- readRDS("dataEdited/binning/depth/bin_depth_clean.rds") %>%
  left_join(MG.metadata %>% select(metagenomeID, date, RM, depth)) %>%
  left_join(metabolic.data)
rm(metabolic.data, MG.metadata)


#### Generate plot for 2017 ####
depth.data %>%
  filter(year(date) == 2017) %>%
  ggplot(aes(y = coverage,
             x = depth,
             group = binID)) +
  geom_point(aes(color = HMS)) +
  geom_line(aes(color = HMS)) +
  scale_color_manual(values = HMS.colors.vector) +
  facet_wrap(~metabolic_assignment + RM, nrow = 4) +
  coord_flip(xlim = c(80, 0)) +
  theme_bw()


#### Generate plot for 2018 ####
depth.data %>%
  filter(year(date) == 2018) %>%
  ggplot(aes(y = coverage,
             x = depth,
             group = binID)) +
  geom_point(aes(color = HMS)) +
  geom_line(aes(color = HMS)) +
  scale_color_manual(values = HMS.colors.vector) +
  facet_wrap(~metabolic_assignment + RM, nrow = 3) +
  coord_flip(xlim = c(80, 0)) +
  theme_bw()


#### Generate plot for 2019 ####
depth.data %>%
  filter(year(date) == 2019) %>%
  filter(RM %in% c(300, 310)) %>%
  ggplot(aes(y = coverage,
             x = depth,
             group = binID)) +
  geom_point(aes(color = HMS)) +
  geom_line(aes(color = HMS)) +
  scale_color_manual(values = HMS.colors.vector) +
  facet_wrap(~metabolic_assignment + RM, nrow = 2) +
  coord_flip(xlim = c(60, 0)) +
  theme_bw()

