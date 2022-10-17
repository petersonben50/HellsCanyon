#### code/metabolic_analyses/nirS_tree.R ####
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


#### Dataset with refs: Read in tree ####
nirS.tree.unrooted <- read.newick(file = "dataEdited/metabolic_analyses/N/nirS/nirS_for_tree.tree")
nirS.tree.rooted <- midpoint(nirS.tree.unrooted)

#### Dataset with refs: Visualize rooted tree ####
pdf("~/Downloads/nirS_tree.pdf",
    height = 80,
    width = 15)
ggtree(nirS.tree.rooted,
       aes(x = 0, xend = 3)) +
  geom_tiplab(size = 2) +
  geom_text2(aes(subset=!isTip, label=node))
dev.off()



#### Final RAxML tree: Read in tree ####
nirS.raxml.unrooted <- read.newick(file = "dataEdited/metabolic_analyses/N/nirS/RAxML_bipartitions.nirS")
# Final RAxML tree: Rename branches
nirS.raxml.unrooted$tip.label[which(nirS.raxml.unrooted$tip.label %in% names(nirS.metadata.vector))] <- nirS.metadata.vector[nirS.raxml.unrooted$tip.label[which(nirS.raxml.unrooted$tip.label %in% names(nirS.metadata.vector))]]

pdf("~/Downloads/nirS_raxml_raw.pdf",
    height = 100,
    width = 15)
ggtree(nirS.raxml.unrooted,
       aes(x = 0, xend = 3)) +
  geom_tiplab(size = 2) +
  geom_text2(aes(subset=!isTip, label=node))
dev.off()


#### Final RAxML tree: Root tree ####
nirS.tree.rooted <- root(nirS.raxml.unrooted,node = 330)


#### Final RAxML tree: Visualize rooted tree ####
pdf("results/metabolic_analyses/nirS_tree.pdf",
    height = 80,
    width = 15)
ggtree(nirS.tree.rooted,
       aes(x = 0, xend = 3)) +
  geom_tiplab(size = 2) + 
  geom_nodelab(aes(x = branch),
               vjust = -.4,
               hjust = 0.6,
               size = 2)
dev.off()
