
#### code/hgcA_analysis/hgcA_data_aggregation.R ####
# Benjamin D. Peterson

rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(readxl)
library(tidyverse)


#### Generate data frame of seqs ####
hgcA.list <- readLines("dataEdited/hgcA_analysis/identification/hgcA_good.txt")
hgcA.df <- data.frame(seqID = hgcA.list,
                      scaffoldID = paste(hgcA.list %>% strsplit("_") %>% sapply("[",1),
                                         hgcA.list %>% strsplit("_") %>% sapply("[",2),
                                         sep = "_"))
rm(hgcA.list)


#### Include truncated info ####
truncated.hgcA.list <- readLines("dataEdited/hgcA_analysis/identification/truncated_hgcA_seq.txt")
hgcA.df <- hgcA.df %>%
  mutate(truncated = (seqID %in% truncated.hgcA.list))
rm(truncated.hgcA.list)


#### Include hgcB info ####
hgcB.list <- readLines("dataEdited/hgcA_analysis/hgcB/hgcB.txt")
hgcB.df <- data.frame(hgcB_ID = hgcB.list,
                      scaffoldID = paste(hgcB.list %>% strsplit("_") %>% sapply("[",1),
                                         hgcB.list %>% strsplit("_") %>% sapply("[",2),
                                         sep = "_"))
hgcB.data.to.read.out <- full_join(hgcA.df,
                                   hgcB.df)
write.csv(hgcB.data.to.read.out,
          "dataEdited/hgcA_analysis/hgcB/hgcB_notes.csv",
          row.names = FALSE)
rm(hgcB.list, hgcB.df, hgcB.data.to.read.out)

hgcB.df <- read_xlsx("dataEdited/hgcA_analysis/hgcB/hgcB_notes.xlsx")
hgcA.df <- hgcA.df %>%
  full_join(hgcB.df)
rm(hgcB.df)


#### Add classification data ####
hgcA.classification <- read.csv("dataEdited/hgcA_analysis/classification/hgcA_taxonomy_table.csv") %>%
  mutate(classification = paste(phylum, class, order,
                                family, genus, species,
                                sep = ";")) %>%
  select(seqID, classification)
hgcA.df <- hgcA.df %>%
  full_join(hgcA.classification)
rm(hgcA.classification)


#### Add cluster IDs ####
clustering.info <- read.table("dataEdited/hgcA_analysis/hgcA_good_acrossYear.tsv",
                              header = TRUE) %>%
  mutate(seqID = id) %>%
  select(seqID, clstr)
hgcA.df <- hgcA.df %>%
  full_join(clustering.info)
rm(clustering.info)



#### Read out data for dereplication ####
write.csv(hgcA.df,
          "dataEdited/hgcA_analysis/hgcA_dereplication_data.csv",
          row.names = FALSE)
rm(hgcA.df)



#### Read data back in ####
hgcA.df <- read_xlsx("dataEdited/hgcA_analysis/hgcA_dereplication.xlsx")

hgcA.df %>%
  filter(representative == TRUE) %>%
  select(seqID) %>%
  unlist(use.names = FALSE) %>%
  writeLines("dataEdited/hgcA_analysis/hgcA_rep_list.txt")

hgcA.df %>%
  filter(usedForAbundance == TRUE) %>%
  select(seqID) %>%
  unlist(use.names = FALSE) %>%
  writeLines("dataEdited/hgcA_analysis/hgcA_repAbundance_list.txt")


#### Add in other info about hgcA ####
hgcA.manual.taxonomy <- read_xlsx("dataEdited/hgcA_analysis/phylogeny/manual_taxonomy.xlsx") %>%
  left_join(hgcA.df %>% select(seqID, clstr)) %>%
  select(clstr, manual_classification)
hgcA.df <- hgcA.df %>%
  left_join(hgcA.manual.taxonomy)


saveRDS(hgcA.df,
        "dataEdited/hgcA_analysis/hgcA_information.rds")
