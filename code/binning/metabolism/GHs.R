#### code/binning/metabolism/GHs.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(readxl)
library(tidyverse)


#### Read in GH data ####
GH.data <- read.table("dataEdited/bins/binAnalysis/metabolism/GHs/cazyme_output.tsv",
                      sep = '\t',
                      skip = 1,
                      col.names = c("geneID", "EC", "GH_ID", "eCAMI", "DIAMOND", "signalP", "NumberOfTools"),
                      header = TRUE) %>%
  select(geneID, GH_ID) %>%
  mutate(GH_class = substr(GH_ID, 1, 2))


#### Read in list of bins to use ####
binList <- readLines("dataEdited/bins/binning/bins_hgcA_keepers/bins_hgcA_keepers_list.txt")


#### Read in taxonomy data ####
bin.data <- read_xlsx("dataEdited/bins/binning/bins_hgcA/bin_dereplication_data_edits.xlsx",
                      sheet = "bin_dereplication_data_trimmed") %>%
  select(binID, HMS, gtdb_tax) %>%
  filter(binID %in% binList)


#### Read in G2B file ####
G2B.df <- read.table("dataEdited/bins/binAnalysis/bin_ORFs/ORFs_G2B.tsv",
                     col.names = c("geneID", "binID")) %>%
  filter(binID %in% binList)



#### Combine data and aggregate counts ####
GH.count.data <- GH.data %>%
  left_join(G2B.df) %>%
  group_by(binID, GH_class) %>%
  summarize(GH_count = n()) %>%
  ungroup() %>%
  left_join(bin.data) %>%
  spread(key = GH_class,
         value = GH_count)


write.csv(GH.count.data %>% select(-c(HMS, gtdb_tax)),
          file = "dataEdited/bins/binAnalysis/metabolism/GHs/clean_cazyme_data.csv",
          row.names = FALSE)
