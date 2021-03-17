#### code/2019_analysis_assembly/sulfur/dsrA/dsr_tree.R ####
# Benjamin D. Peterson


#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(ape)
library(ggtree)
library(phangorn)
library(readxl)
library(tidyverse)
library(treeio)


#### Read in tree ####
dsrA.tree.unrooted <- read.newick(file = "dataEdited/metabolic_analyses/sulfur/dsrA_phylogeny_masked.tree")
dsrA.tree.unrooted$tip.label[grep("ˆ", dsrA.tree.unrooted$tip.label)] <- gsub("ˆ",
                                                                              "u",
                                                                              dsrA.tree.unrooted$tip.label[grep("ˆ", dsrA.tree.unrooted$tip.label)])


#### Visualize unrooted tree ####
pdf("dataEdited/metabolic_analyses/sulfur/dsrA_tree_original_unrooted.pdf",
    width = 12,
    height = 120)
ggtree(dsrA.tree.unrooted) +
  geom_tiplab(size = 2) +
  geom_text2(aes(subset=!isTip, label=node))
dev.off()


#### Root tree ####
dsrA.tree <- root(dsrA.tree.unrooted,
                  node = 1125)


#### Visualize rooted tree ####
pdf("dataEdited/metabolic_analyses/sulfur/dsrA_tree_original_rooted.pdf",
    width = 12,
    height = 120)
ggtree(dsrA.tree) +
  geom_tiplab(size = 2) +
  geom_text2(aes(subset=!isTip, label=node))
dev.off()

#### Get list of dsrA genes ####
dsrA.list <- readLines("dataEdited/metabolic_analyses/sulfur/dsrA_list.txt")


#### Get list of rdsr genes ####
# First isolate those tip labels
endpoints <- c("Nitrospirae_bacterium_RBG_19FT_COMBO_42_15",
               "RBG_16_scaffold_151951")
mrca.of.interest <- getMRCA(dsrA.tree,
                            endpoints)
rdsr.tree <- tree_subset(dsrA.tree,
                         mrca.of.interest,
                         levels_back = 0)
ggtree(rdsr.tree) +
  geom_tiplab(size = 2)

rdsr.list <- rdsr.tree$tip.label[rdsr.tree$tip.label %in% dsrA.list]
writeLines(rdsr.list,
           "dataEdited/metabolic_analyses/sulfur/dsrA_rev_list.txt")


#### Find reductive dsrA genes ####
red.dsrA <- dsrA.list[which(!(dsrA.list %in% rdsr.list))]


#### Read out list of reductive dsrA genes ####

writeLines(red.dsrA,
           "dataEdited/metabolic_analyses/sulfur/dsrA_red_list.txt")


