
#### Clean up crew on line 5 ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(tidyverse)


#### Read in normalization vector ####
normalized.coverage.vector <- readRDS("dataEdited/scg_abundance/scg_normalization_vector.rds")


#### Read in metadata ####
metadata.df <- read.csv("metadata/metagenome_metadata.csv")


#### Read in depth data ####
list.o.depths <- list.files(path = "dataEdited/bins/binAnalysis/hqBins/depth",
                            pattern = "depth.tsv",
                            full.names = TRUE)

raw.depth.counts <- lapply(list.o.depths,
                           function(fileName) {
                             
                             metagenomeID = fileName %>%
                               gsub("dataEdited/bins/binAnalysis/hqBins/depth/", "", .) %>%
                               gsub("_depth.tsv", "", .)
                             
                             read.csv(fileName,
                                      stringsAsFactors = FALSE) %>%
                               mutate(meanCoverageBin = round(meanCoverageBin * normalized.coverage.vector[metagenomeID], 4)) %>%
                               mutate(read.origin = metagenomeID)
                             
                           })


#### Combine data into dataframe ####
depth.df <- do.call(rbind,
                    raw.depth.counts) %>%
  rename(metagenomeID = read.origin,
         coverage = meanCoverageBin)
#### Write out file
saveRDS(depth.df,
        "dataEdited/bins/binAnalysis/hqBins/bin_depth_clean.rds")
