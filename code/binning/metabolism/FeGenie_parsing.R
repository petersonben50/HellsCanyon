#### code/binning/metabolism/FeGenie_parsing.R ####
# Benjamin D. Peterson


#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(tidyverse)


#### Read in bin data ####
bin.data <- read.csv("dataEdited/binning/bin_data.csv",
                     stringsAsFactors = FALSE) %>%
  select(binID, HMS, assemblyID, gtdb_tax,
         checkM_completeness, checkM_contamination)


#### Read in FeGenie data ####
FeGenie.data <- read.csv('dataEdited/binning/metabolism/FeGenie/All_EET_Proteins_bothMethods_HellsCanyon13genomes.csv') %>%
  rename(binID = MAG)



#### Combine data ####
all.data <- full_join(bin.data,
                      FeGenie.data) %>%
  filter(orf != "Not found")


#### Check out Bacteroidetes ####

bact.data <- all.data %>%
  filter(binID == "anvio_hgcA_0130")
