#### code/phylogenies/PVC_tree.R ####
# Benjamin D. Peterson

# This script will generate images of
# the phylogenetic trees for the PVC bins.



#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon/")
library(ggtree)
library(phangorn)
library(readxl)
library(tidyverse)
library(treeio)
cb.translator <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")

#### Initial check of tree ####

# Read in tree
tree.name <- "dataEdited/binning/phylogeny/PVC/RAxML_bipartitions.PVC_rp16"
hgcA.tree.unrooted <- read.newick(tree.name)
rm(tree.name)

# Check out unrooted tree
pdf("dataEdited/binning/phylogeny/PVC/PVC_tree_FastTree_unrooted.pdf",
    height = 10,
    width = 5)
ggtree(hgcA.tree.unrooted,
       aes(x = 0,
           xend = 2)) + 
  geom_tiplab(size=2.5, align = TRUE)
dev.off()

hgcA.tree <- midpoint(hgcA.tree.unrooted)
ggtree(hgcA.tree,
       aes(x = 0,
           xend = 2)) + 
  geom_tiplab(size=2.5, align = TRUE)
