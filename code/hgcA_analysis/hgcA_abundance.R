#### code/hgcA_analysis/hgcA_abundance.R ####
# Written by Benjamin D. Peterson


#### Get ready. Get set ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(lubridate)
library(patchwork)
library(readxl)
library(tidyverse)


#### Read in hgcA seq ID list ####
hgcA.list <- readLines("dataEdited/hgcA_analysis/hgcA_repAbundance_list.txt")


#### Read in hgcA classification ####
tax.data <- readRDS("dataEdited/hgcA_analysis/hgcA_information.rds") %>%
  select(seqID, manual_classification)


#### Read in data ####
hgcA.data <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv") %>%
  filter(seqID %in% hgcA.list) %>%
  left_join(tax.data)


#### Function to generate plots ####
RM.of.interest <- "RM300"
year.of.interest <- 2018
hgcA.profiling <- function(RM.of.interest,
                           year.of.interest) {
  hgcA.data %>%
    filter(year(date) == year.of.interest,
           RM == RM.of.interest) %>%
    ggplot(aes(x = depth,
               y = coverage,
               fill = manual_classification)) +
    geom_bar(stat = "identity") +
    xlim(80, 0) +
    coord_flip() +
    theme_classic() +
    labs(title = paste("hgcA coverage at ",
                       RM.of.interest,
                       " in ",
                       year.of.interest,
                       sep = ""))
  
}

profile.286.2018 <- hgcA.profiling(RM.of.interest = "RM286",
                                   year.of.interest = 2018)
profile.300.2018 <- hgcA.profiling(RM.of.interest = "RM300",
                                   year.of.interest = 2018)

profile.286.2018 + profile.300.2018



profile.286.2017 <- hgcA.profiling(RM.of.interest = "RM286",
                                   year.of.interest = 2017)
profile.300.2017 <- hgcA.profiling(RM.of.interest = "RM300",
                                   year.of.interest = 2017)

profile.286.2017 + profile.300.2017
