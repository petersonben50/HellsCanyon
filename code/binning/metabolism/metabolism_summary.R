#### code/binning/metabolism/metabolism_summary.R ####
# Benjamin D. Peterson


#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(tidyverse)


#### Read in HMM data ####
metabolism.data <- read.table("dataEdited/binning/metabolism/all_bin_counts.tsv",
                              header = TRUE) %>%
  spread(key = proteinName,
         value = counts)


#### Read in bin data ####
bin.data <- read.csv("dataEdited/binning/bin_data.csv",
                     stringsAsFactors = FALSE) %>%
  select(binID, HMS, assemblyID, gtdb_tax,
         checkM_completeness, checkM_contamination)


#### Read in MHC data ####
MHC.data <- read.table('dataEdited/binning/metabolism/MHCs/heme_count_bins.tsv',
                       sep = '\t',
                       header = TRUE) %>%
  group_by(binID) %>%
  summarise(MHC_count = n())


#### Join data together ####
all.data <- left_join(bin.data,
                      MHC.data) %>%
  left_join(metabolism.data)
# Add EET info
all.data <- all.data %>%
  mutate(PCC = (binID == "anvio_hgcA_0210")*1)


#### Save out data ####
write.csv(all.data,
          "dataEdited/binning/metabolism/metabolic_summary.csv",
          row.names = FALSE)
