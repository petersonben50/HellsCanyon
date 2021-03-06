#### code/hgcA_analysis/clean_hgcA_coverage.R ####
# Written by Benjamin D. Peterson


#### Get ready. Get set ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(readxl)
library(stringr)
library(tidyverse)



#### Read in normalization vector ####
# Normalized by SCG coverage
normalized.coverage.vector <- readRDS("dataEdited/scg_abundance/scg_normalization_vector.rds")
# normalized to 10 million reads
# normalized.coverage.vector <- readRDS("dataEdited/metagenomes/reports/metagenome_normalization_vector.rds")


#### List of file names with depth info ####
list.o.depths <- list.files(path = "dataEdited/hgcA_analysis/depth",
                            pattern = "hgcA_depth",
                            full.names = TRUE)


#### Metadata vector ####
metadata.df <- read.csv("metadata/metagenome_metadata.csv")


#### Read in and normalize depth data ####

raw.depth.counts <- lapply(list.o.depths,
                           function(fileName) {
                             
                             metagenomeID = fileName %>%
                               gsub("dataEdited/hgcA_analysis/depth/", "", .) %>%
                               gsub("_hgcA_depth.tsv", "", .)
                             
                             read.table(fileName,
                                        stringsAsFactors = FALSE,
                                        col.names = c("seqID", "depth", "length")) %>%
                               mutate(depth = depth * normalized.coverage.vector[metagenomeID]) %>%
                               mutate(read.origin = metagenomeID)
                             
                             })
rm(list.o.depths,
   normalized.coverage.vector)



#### Combine data into dataframe ####
coverage.all.sequences <- do.call(rbind,
                                  raw.depth.counts) %>%
  spread(key = read.origin, value = depth, fill = NA) %>%
  select(-length) %>%
  gather(key = "metagenomeID",
         value = "coverage",
         -1) %>%
  filter(!is.na(coverage))

#### Add metadata ####
coverage.all.sequences <- coverage.all.sequences %>%
  left_join(metadata.df) %>%
  rename(scaffoldID = seqID)
rm(raw.depth.counts,
   metadata.df)


#### Add hgcA seq ID ####
hgcA.df <- read_xlsx("dataEdited/hgcA_analysis/hgcA_dereplication.xlsx") %>%
  select(seqID, scaffoldID)
coverage.all.sequences <- full_join(hgcA.df,
                                    coverage.all.sequences)

#### Read out dataframe ####
write.csv(coverage.all.sequences,
          "dataEdited/hgcA_analysis/depth/hgcA_coverage.csv",
          row.names = FALSE)
