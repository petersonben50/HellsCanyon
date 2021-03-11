#### code/hgcA_analysis/hgcA_abundance.R ####
# Written by Benjamin D. Peterson


#### Get ready. Get set ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(lubridate)
library(patchwork)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
cb.translator <- c(cb.translator, "gray50")
names(cb.translator)[length(cb.translator)] <- "gray50"

#### Read in hgcA seq ID list ####
hgcA.list <- readLines("dataEdited/hgcA_analysis/hgcA_repAbundance_list.txt")


#### Read in hgcA classification ####
tax.data <- readRDS("dataEdited/hgcA_analysis/hgcA_information.rds") %>%
  select(seqID, manual_classification)


#### Make color vector ####
hgcA.manual.taxonomy <- read_xlsx("dataEdited/hgcA_analysis/phylogeny/manual_taxonomy.xlsx",
                                  sheet = "colors_to_use")
color.vector <- cb.translator[hgcA.manual.taxonomy$colorsToUse]
names(color.vector) <- hgcA.manual.taxonomy$seqID

#### Read in data ####
hgcA.data <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv") %>%
  filter(seqID %in% hgcA.list) %>%
  left_join(tax.data)


#### Function to generate plots ####
RM.of.interest <- "RM300"
year.of.interest <- 2018
hgcA.profiling <- function(data.to.use = hgcA.data,
                           RM.of.interest,
                           year.of.interest,
                           legend.position.to.use = "default",
                           xlim.to.use = c(80, 0)) {
  graph.to.make <- data.to.use %>%
    filter(year(date) == year.of.interest,
           RM == RM.of.interest) %>%
    ggplot(aes(x = depth,
               y = coverage,
               fill = manual_classification)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color.vector) +
    xlim(xlim.to.use) +
    coord_flip() +
    theme_classic() +
    labs(title = paste("hgcA coverage at ",
                       RM.of.interest,
                       " in ",
                       year.of.interest,
                       sep = ""))
  
  if (legend.position.to.use[1] != "default") {
    graph.to.make <- graph.to.make +
      theme(legend.position = legend.position.to.use)
  }
  
  graph.to.make
  
}

profile.286.2018 <- hgcA.profiling(RM.of.interest = "RM286",
                                   year.of.interest = 2018,
                                   legend.position.to.use = "none")
profile.300.2018 <- hgcA.profiling(RM.of.interest = "RM300",
                                   year.of.interest = 2018,
                                   legend.position.to.use = c(0.7, 0.8))

profile.286.2018 + profile.300.2018



profile.286.2017 <- hgcA.profiling(RM.of.interest = "RM286",
                                   year.of.interest = 2017,
                                   legend.position.to.use = "none")
profile.300.2017 <- hgcA.profiling(RM.of.interest = "RM300",
                                   year.of.interest = 2017,
                                   legend.position.to.use = c(0.7, 0.8))

profile.286.2017 + profile.300.2017


# Zoom in on metalimnion in 2017
profile.286.2017 <- hgcA.profiling(data.to.use = hgcA.data %>%
                                     filter(depth < 41),
                                   RM.of.interest = "RM286",
                                   year.of.interest = 2017,
                                   legend.position.to.use = c(0.7, 0.8),
                                   xlim.to.use = c(41, 0))
profile.286.2017