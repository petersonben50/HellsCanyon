#### code/binning/PVC_details/PVC_hgcA_alignment.R ####
# Benjamin D. Peterson

# This script will generate images of
# the phylogenetic trees for the PVC bins.


#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon/")
library(Biostrings)
library(DECIPHER)


#### Investigate alignment ####
readAAStringSet("dataEdited/bins/binAnalysis/PVC_details/hgcA_analysis/hgcA_hits/PVC_hgcA.afa",
                format="fasta",
                nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE) %>%
  BrowseSeqs(htmlFile = "dataEdited/bins/binAnalysis/PVC_details/hgcA_analysis/hgcA_hits/PVC_hgcA.html",
             colorPatterns = c(226, 234))
