#### code/binning/phylogenies/prolixibacter_tree.R ####
# Benjamin D. Peterson

# This script will generate images of
# the phylogenetic trees for the bacteroidetes bins.


#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon/")
library(ggtree)
library(phangorn)
library(readxl)
library(tidyverse)
library(treeio)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Investigate GTDB tree ####
GTDB.tree <- read.newick("dataEdited/bins/binAnalysis/phylogeny/gtdbtk.bac120.classify.tree")

# Subset tree
prolix.tree <- tree_subset(GTDB.tree,
                         "anvio_hgcA_0130",
                         levels_back = 2)
ggtree(prolix.tree,
       aes(x = 0,
           xend = 0.4)) +
  geom_tiplab(size = 2)
prolix.tree$tip.label
