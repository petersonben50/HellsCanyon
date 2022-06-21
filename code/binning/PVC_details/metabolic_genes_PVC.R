#### code/binning/PVC_details/metabolic_genes_PVC.R ####
# Benjamin D. Peterson


#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon/")
library(readxl)
library(tidyverse)
library(treeio)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


####-------------Read in data-------------####
pvc.tree <- read.newick("dataEdited/bins/binAnalysis/PVC_details/tree_generation/RAxML_bipartitions.HCC_PVC")
bin.list <- pvc.tree$tip.label
rm(pvc.tree)

MHC.data <- read.table("dataEdited/bins/binAnalysis/PVC_details/metabolism/heme_count_bins.tsv",
                       sep = '\t', header = TRUE) %>%
  group_by(binID) %>%
  summarise(MHC_count = n())
HMM.data <- read.table("dataEdited/bins/binAnalysis/PVC_details/metabolism/all_bin_counts.tsv",
                       sep = '\t', header = TRUE) %>%
  spread(key = proteinName,
         value = counts,
         fill = 0)
all.data <- full_join(MHC.data, HMM.data) %>%
  filter(binID %in% bin.list)


####-------------Read in data-------------####
MHC.cutoff <- 5
high.mhc.bins <- all.data %>%
  filter(MHC_count > MHC.cutoff) %>%
  select(binID) %>%
  unlist(use.names = FALSE)
saveRDS(high.mhc.bins, file = "dataEdited/bins/binAnalysis/PVC_details/metabolism/high_MHC_bins.rds")

low.mhc.bins <- all.data$binID[!(all.data$binID %in% high.mhc.bins)]
saveRDS(low.mhc.bins, file = "dataEdited/bins/binAnalysis/PVC_details/metabolism/low_MHC_bins.rds")
