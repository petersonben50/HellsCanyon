#### code/binning/PVC_details/hgcB_alignment.R ####
# Benjamin D. Peterson

# This script will generate images of
# the phylogenetic trees for the PVC bins.


#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon/")
library(Biostrings)
library(DECIPHER)
library(tidyverse)


#### Investigate alignment for hgcB hits ####
readAAStringSet("dataEdited/bins/binAnalysis/PVC_details/hgcA_analysis/hgcB/PVC_hgcB.afa",
                format="fasta",
                nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE) %>%
  BrowseSeqs(htmlFile = "dataEdited/bins/binAnalysis/PVC_details/hgcA_analysis/hgcB/hgcB_hits.html",
             colorPatterns = c(51,57))


#### Investigate alignment for hgcB hits ####
readAAStringSet("dataEdited/bins/binAnalysis/PVC_details/hgcA_analysis/hgcB/hgcB_suspects_interrogation.afa",
                format="fasta",
                nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE) %>%
  BrowseSeqs(htmlFile = "dataEdited/bins/binAnalysis/PVC_details/hgcA_analysis/hgcB/hgcB_suspects_interrogation.html",
             colorPatterns = c(61,66))
