#### code/metabolic_analyses/narG_tree.R ####
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


#### References: Read in tree ####
narG.tree.unrooted <- read.newick(file = "references/narG/narG_luke_database_fastTree.tree")
narG.tree.rooted <- root(narG.tree.unrooted,node = 435)


#### References: Read in metadata ####
narG.metadata <- read_xlsx("references/narG/narG_luke_database.xlsx",
                           sheet = "narG_luke_metadata")
narG.metadata.vector <- paste(narG.metadata$treeID,
                              narG.metadata$description,
                              sep = " - ")
names(narG.metadata.vector) <- narG.metadata$treeID

#### References: Rename branches ####
narG.tree.rooted$tip.label[which(narG.tree.rooted$tip.label %in% names(narG.metadata.vector))] <- narG.metadata.vector[narG.tree.rooted$tip.label[which(narG.tree.rooted$tip.label %in% names(narG.metadata.vector))]]

#### References: Visualize rooted tree ####
ggtree(narG.tree.rooted,
       aes(x = 0, xend = 3)) +
  geom_tiplab(size = 2) +
  geom_text2(aes(subset=!isTip, label=node))
  


#### Dataset with refs: Read in tree ####
narG.tree.unrooted <- read.newick(file = "dataEdited/metabolic_analyses/N/narG/narG_for_tree.tree")
narG.tree.rooted <- root(narG.tree.unrooted,node = 435)


#### Dataset with refs: Read in metadata ####
narG.metadata <- read_xlsx("references/narG/narG_luke_database.xlsx",
                           sheet = "narG_luke_metadata")
narG.metadata.vector <- paste(narG.metadata$treeID,
                              narG.metadata$description,
                              sep = " - ")
names(narG.metadata.vector) <- narG.metadata$treeID

#### Dataset with refs: Rename branches ####
narG.tree.rooted$tip.label[which(narG.tree.rooted$tip.label %in% names(narG.metadata.vector))] <- narG.metadata.vector[narG.tree.rooted$tip.label[which(narG.tree.rooted$tip.label %in% names(narG.metadata.vector))]]

#### Dataset with refs: Visualize rooted tree ####
ggtree(narG.tree.rooted,
       aes(x = 0, xend = 3)) +
  geom_tiplab(size = 2) +
  geom_text2(aes(subset=!isTip, label=node))




#### Final RAxML tree: Read in tree ####
narG.raxml.unrooted <- read.newick(file = "dataEdited/metabolic_analyses/N/narG/RAxML_bipartitions.narG")
# Final RAxML tree: Rename branches
narG.raxml.unrooted$tip.label[which(narG.raxml.unrooted$tip.label %in% names(narG.metadata.vector))] <- narG.metadata.vector[narG.raxml.unrooted$tip.label[which(narG.raxml.unrooted$tip.label %in% names(narG.metadata.vector))]]

pdf("~/Downloads/narG_raxml_raw.pdf",
    height = 100,
    width = 15)
ggtree(narG.raxml.unrooted,
       aes(x = 0, xend = 3)) +
  geom_tiplab(size = 2) +
  geom_text2(aes(subset=!isTip, label=node))
dev.off()


#### Final RAxML tree: Root tree ####
narG.tree.rooted <- root(narG.raxml.unrooted,node = 330)


#### Final RAxML tree: Visualize rooted tree ####
pdf("results/metabolic_analyses/narG_tree.pdf",
    height = 80,
    width = 15)
ggtree(narG.tree.rooted,
       aes(x = 0, xend = 3)) +
  geom_tiplab(size = 2) + 
  geom_nodelab(aes(x = branch),
               vjust = -.4,
               hjust = 0.6,
               size = 2)
dev.off()
