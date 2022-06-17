#### code/binning/metabolism/metabolism_summary.R ####
# Benjamin D. Peterson


#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(readxl)
library(tidyverse)


#### Read in HMM data ####
metabolism.data <- read.table("dataEdited/bins/binAnalysis/metabolism/batch_HMMS/all_bin_counts.tsv",
                              header = TRUE) %>%
  spread(key = proteinName,
         value = counts)
write.csv(metabolism.data,
          "dataEdited/bins/binAnalysis/metabolism/batch_HMMS/all_bin_counts_spread.csv",
          row.names = FALSE)

#### Read in bin data ####
bin.data <- read_xlsx("dataEdited/bins/binning/bins_hgcA/bin_dereplication_data_edits.xlsx",
                      sheet = "bin_dereplication_data_trimmed") %>%
  select(binID, HMS, assemblyID, gtdb_tax,
         checkM_completeness, checkM_contamination)


#### Read in MHC data ####
MHC.data <- read.table('dataEdited/bins/binAnalysis/metabolism/MHCs/heme_count_bins.tsv',
                       sep = '\t',
                       header = TRUE) %>%
  group_by(binID) %>%
  summarise(MHC_count = n())
write.csv(MHC.data,
          'dataEdited/bins/binAnalysis/metabolism/MHCs/MHC_count_bins.csv',
          row.names = FALSE)

#### Join data together ####
all.data <- left_join(bin.data,
                      MHC.data) %>%
  left_join(metabolism.data)


# Add EET info
# all.data <- all.data %>%
#   mutate(PCC = (binID == "anvio_hgcA_0210")*1)


#### Save out data ####
write.csv(all.data,
          "dataEdited/binning/metabolism/metabolic_summary.csv",
          row.names = FALSE)
