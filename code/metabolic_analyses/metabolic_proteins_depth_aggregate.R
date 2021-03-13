#### code/metabolic_analyses/metabolic_proteins_depth_aggregate.R ####
# Benjamin D. Peterson


#### Clean up crew on line 5 ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(tidyverse)


#### Read in normalization vector ####
normalized.coverage.vector <- readRDS("dataEdited/scg_abundance/scg_normalization_vector.rds")


#### Read in metadata ####
metadata.df <- read.csv("metadata/metagenome_metadata.csv")


#### Read in metabolic gene key ####
metabolic.gene.key <- read.csv("dataEdited/metabolic_analyses/metabolic_gene_key.csv") %>%
  mutate(scaffoldID = paste(geneID %>% strsplit("_") %>% sapply("[", 1),
                            geneID %>% strsplit("_") %>% sapply("[", 2),
                            sep = "_")) %>%
  select(scaffoldID, geneName)


#### Read in depth data ####
list.o.depths <- list.files(path = "dataEdited/metabolic_analyses/depth",
                            pattern = "depth.tsv",
                            full.names = TRUE)

raw.depth.counts <- lapply(list.o.depths,
                           function(fileName) {
                             
                             metagenomeID = fileName %>%
                               gsub("dataEdited/metabolic_analyses/depth/", "", .) %>%
                               gsub("_depth.tsv", "", .)
                             
                             read.table(fileName,
                                        stringsAsFactors = FALSE,
                                        header = FALSE,
                                        col.names = c("scaffoldID", "depth")) %>%
                               mutate(depth = round(depth * normalized.coverage.vector[metagenomeID], 4)) %>%
                               mutate(read.origin = metagenomeID)
                             
                           })


#### Combine data into dataframe ####
depth.df <- do.call(rbind,
                    raw.depth.counts) %>%
  rename(metagenomeID = read.origin,
         coverage = depth)


#### Add in metadata ####
all.data <- metabolic.gene.key %>%
  full_join(depth.df) %>%
  left_join(metadata.df)


#### Write out file
saveRDS(all.data,
        "dataEdited/metabolic_analyses/depth/metabolicProtein_depth_clean.rds")
