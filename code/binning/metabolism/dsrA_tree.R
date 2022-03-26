#### code/binning/metabolism/dsr_tree.R ####
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
dsrA.tree.unrooted <- read.newick(file = "dataEdited/bins/binAnalysis/metabolism/batch_HMMS/dsrA/dsrA_phylogeny_trimmed.tree")
dsrA.tree.unrooted$tip.label[grep("ˆ", dsrA.tree.unrooted$tip.label)] <- gsub("ˆ",
                                                                              "u",
                                                                              dsrA.tree.unrooted$tip.label[grep("ˆ", dsrA.tree.unrooted$tip.label)])


#### Visualize unrooted tree ####
pdf("dataEdited/bins/binAnalysis/metabolism/batch_HMMS/dsrA/dsrA_tree_original_unrooted.pdf",
    width = 12,
    height = 80)
ggtree(dsrA.tree.unrooted) +
  geom_tiplab(size = 2) +
  geom_text2(aes(subset=!isTip, label=node))
dev.off()


#### Root tree ####
dsrA.tree <- root(dsrA.tree.unrooted,
                  node = 771)


#### Get list of dsrA genes, color bin genes in red ####
dsrA.list <- readLines("dataEdited/bins/binAnalysis/metabolism/batch_HMMS/dsrA/dsrA_list.txt") %>%
  strsplit(" ") %>% sapply("[", 1)
tip.label.color <- rep("black",
                       length(dsrA.tree$tip.label))
tip.label.color[dsrA.tree$tip.label %in% dsrA.list] <- 'red'

#### Visualize rooted tree ####
pdf("dataEdited/bins/binAnalysis/metabolism/batch_HMMS/dsrA/dsrA_tree_original.pdf",
    width = 12,
    height = 120)
ggtree(dsrA.tree) +
  geom_tiplab(size = 2,
              colour = tip.label.color) +
  geom_text2(aes(subset=!isTip, label = node))
dev.off()


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


#### Find reverse dsrA genes ####
rdsr.tree$tip.label[rdsr.tree$tip.label %in% dsrA.list]
# They're all reductive
