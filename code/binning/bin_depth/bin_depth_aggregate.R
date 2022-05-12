#### code/metabolic_analyses/bin_depth_aggregate.R ####
# Benjamin D. Peterson


#### Clean up crew on line 5 ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(readxl)
library(tidyverse)


#### Read in normalization vector ####
normalized.coverage.vector <- readRDS("dataEdited/scg_abundance/scg_normalization_vector.rds")


#### Read in metadata ####
metadata.df <- read.csv("metadata/metagenome_metadata.csv")


#### Read in depth data ####
list.o.depths <- list.files(path = "dataEdited/bins/binAnalysis/depth",
                            pattern = "depth.tsv",
                            full.names = TRUE)

raw.depth.counts <- lapply(list.o.depths,
                           function(fileName) {
                             
                             metagenomeID = fileName %>%
                               gsub("dataEdited/bins/binAnalysis/depth/", "", .) %>%
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


# See what percent of microbial community from assembly was
# binned as hgcA+ bins

#### Read in taxonomy data ####
bin.data <- read_xlsx("dataEdited/bins/binning/bins_hgcA/bin_dereplication_data_edits.xlsx",
                      sheet = "bin_dereplication_data_trimmed") %>%
  select(binID, HMS, gtdb_tax)

depth.df.tax <- depth.df %>%
  left_join(bin.data) %>%
  group_by(metagenomeID, HMS, gtdb_tax) %>%
  summarize(coverage = mean(coverage)) %>%
  left_join(metadata.df)

test <- depth.df.tax %>%
  group_by(metagenomeID) %>%
  summarise(total_coverage = sum(coverage)) %>%
  left_join(metadata.df)


#### Write out file
saveRDS(depth.df,
        "dataEdited/bins/binAnalysis/depth/bin_depth_clean.rds")
