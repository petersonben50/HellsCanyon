#### code/binning/geobacterales/geobacterales_derep.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(tidyverse)


#### Read in taxonomy data ####
geo.tax.data <- read_xlsx("dataEdited/bins/binning/autoBinning/taxonomy_summary.xlsx") %>%
  mutate(family = classification %>%
           strsplit(";f__") %>% sapply("[", 2) %>%
           strsplit(";g__") %>% sapply("[", 1),
         genus = classification %>%
           strsplit(";g__") %>% sapply("[", 2) %>%
           strsplit(";s__") %>% sapply("[", 1)) %>%
  filter(family == "Pelobacteraceae")
