#### code/hgcA_analysis/clean_hgcA_tree.R ####
# Benjamin D. Peterson

# This script will look at our phylogenetic tree 
# iterations and ultimately clean up the final one.


#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon/")
library(ggtree)
library(phangorn)
library(readxl)
library(tidyverse)
library(treeio)
cb.translator <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")




#### Read in tree ####
hgcA.tree <- readRDS("dataEdited/hgcA_analysis/phylogeny/hgcA_clean_tree.rds")


#### Save out color vector ####
color.vector <- readRDS("dataEdited/hgcA_analysis/phylogeny/hgcA_clean_tree_color_vector.rds")


#### Rename one tip ####
hgcA.tree$tip.label[grep("fall2017coassembly_000000072407_5",
                         hgcA.tree$tip.label)] <- paste(hgcA.tree$tip.label[grep("fall2017coassembly_000000072407_5",
                                                                                 hgcA.tree$tip.label)],
                                                        " - anvio_hgcA_0070",
                                                        sep = "")


#### Remove bootstraps < 50 ####
hgcA.tree$node.label[hgcA.tree$node.label == "Root"] <- ""
hgcA.tree$node.label[as.numeric(hgcA.tree$node.label) < 50] <- ""


#### Visualize tree ####
pdf("results/manuscript_figures/geochem_figs_for_supp_text/hgcA_tree.pdf",
    height = 60,
    width = 15)
ggtree(hgcA.tree, aes(x = 0,
                      xend = 9)) + 
  geom_tiplab(size=2.5, align = TRUE, col = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 2) +
  geom_treescale(x = 0.05,
                 y = 30,
                 width = 1)
dev.off()

