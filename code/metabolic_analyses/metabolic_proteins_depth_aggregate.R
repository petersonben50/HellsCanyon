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


#### Add elevation data ####
metadata.df <- metadata.df %>%
  left_join(read.csv("dataEdited/geochem/geochem_WC.csv") %>%
              select(RM, date, depth, elevation_m) %>%
              unique())


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


#### Adjust dsrA labels ####
rdsrA.list <- readLines("dataEdited/metabolic_analyses/sulfur/dsrA_rev_list.txt")
rdsrA.scaffolds <- paste(rdsrA.list %>% strsplit("_") %>% sapply("[", 1),
                         rdsrA.list %>% strsplit("_") %>% sapply("[", 2),
                         sep = "_")
dsrA.list <- readLines("dataEdited/metabolic_analyses/sulfur/dsrA_red_list.txt")
dsrA.scaffolds <- paste(dsrA.list %>% strsplit("_") %>% sapply("[", 1),
                        dsrA.list %>% strsplit("_") %>% sapply("[", 2),
                        sep = "_")
all.data[which((all.data$scaffoldID %in% rdsrA.scaffolds) & (all.data$geneName == 'dsrA')), "geneName"] <- "rdsrA"
all.data[which((all.data$scaffoldID %in% dsrA.scaffolds) & (all.data$geneName == 'dsrA')), "geneName"]



#### Add EET gene data ####
EET.data <- readRDS("dataEdited/metabolic_analyses/BBOMP/bbomp_depth_clean.rds") %>%
  left_join(read.csv("dataEdited/geochem/geochem_WC.csv") %>%
              select(RM, date, depth, elevation_m) %>%
              unique())
all.data.w.EET <- rbind(all.data,
                        EET.data) %>%
  filter(!is.na(geneName))


#### Write out file ####
saveRDS(all.data.w.EET,
        "dataEdited/metabolic_analyses/depth/metabolicProtein_depth_clean.rds")



#### Save out summarized file ####
summarized.gene.data <- all.data.w.EET %>%
  group_by(metagenomeID, date, RM, depth, elevation_m, geneName) %>%
  summarize(coverage = sum(coverage)) %>%
  spread(key = geneName,
         value = coverage)
write.csv(x = summarized.gene.data,
          file = "dataEdited/metabolic_analyses/summarized_gene_data.csv",
          row.names = FALSE)

