#### code/mer/merB_abundance.R ####
# Written by Benjamin D. Peterson


#### Get ready. Get set ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(lubridate)
library(patchwork)
library(readxl)
library(tidyverse)


#### Read in merB seq ID list ####
merB.list <- readLines("dataEdited/mer/dereplication/merB_derep_list.txt")
merB.df <- data.frame(seqID = merB.list,
                      scaffoldID = paste(merB.list %>% strsplit("_") %>% sapply("[", 1),
                                         merB.list %>% strsplit("_") %>% sapply("[", 2),
                                         sep = "_"))

#### Read in merB classification ####
class.data <- read.csv("dataEdited/mer/identification/final_groups/group_key.csv")
merB.df <- merB.df %>%
  inner_join(class.data)


#### Read in data ####
merB.data <- read.csv("dataEdited/mer/depth/merB_coverage.csv") %>%
  left_join(merB.df) %>%
  filter(seqID %in% merB.list)


#### Function to generate plots ####
RM.of.interest <- 300
year.of.interest <- 2018
merB.profiling <- function(RM.of.interest,
                           year.of.interest,
                           cluster.of.interest = NULL,
                           where.to.put.legend = "none") {
  if (!is.null(cluster.of.interest)) {
    merB.data <- merB.data %>%
      filter(clusterName == cluster.of.interest)
  }
  
  merB.data %>%
    filter(year(date) == year.of.interest,
           RM == RM.of.interest) %>%
    ggplot(aes(x = depth,
               y = coverage,
               fill = clusterName)) +
    geom_bar(stat = "identity") +
    xlim(80, 0) +
    coord_flip() +
    theme_classic() +
    labs(title = paste("merB coverage at ",
                       RM.of.interest,
                       " in ",
                       year.of.interest,
                       sep = "")) +
    theme(legend.position = where.to.put.legend)
  
}



profile.286.2017 <- merB.profiling(RM.of.interest = 286,
                                   year.of.interest = 2017)
profile.300.2017 <- merB.profiling(RM.of.interest = 300,
                                   year.of.interest = 2017)

profile.286.2017 + profile.300.2017

profile.286.2018 <- merB.profiling(RM.of.interest = 286,
                                   year.of.interest = 2018)
profile.300.2018 <- merB.profiling(RM.of.interest = 300,
                                   year.of.interest = 2018)

profile.286.2018 + profile.300.2018



profile.300.2019 <- merB.profiling(RM.of.interest = 300,
                                   year.of.interest = 2019)
profile.310.2019 <- merB.profiling(RM.of.interest = 310,
                                   year.of.interest = 2019)
profile.300.2019 + profile.310.2019



profile.286.2018 <- merB.profiling(RM.of.interest = 286,
                                   year.of.interest = 2018,
                                   cluster.of.interest = "merB3")
profile.300.2018 <- merB.profiling(RM.of.interest = 300,
                                   year.of.interest = 2018,
                                   cluster.of.interest = "merB3")

profile.286.2018 + profile.300.2018



profile.286.2017 <- merB.profiling(RM.of.interest = 286,
                                   year.of.interest = 2017,
                                   cluster.of.interest = "merB3")
profile.300.2017 <- merB.profiling(RM.of.interest = 300,
                                   year.of.interest = 2017,
                                   cluster.of.interest = "merB3")

profile.286.2017 + profile.300.2017



profile.286.2017 <- merB.profiling(RM.of.interest = 286,
                                   year.of.interest = 2017,
                                   cluster.of.interest = "merB_fused")
profile.300.2017 <- merB.profiling(RM.of.interest = 300,
                                   year.of.interest = 2017,
                                   cluster.of.interest = "merB_fused")

profile.286.2017 + profile.300.2017
