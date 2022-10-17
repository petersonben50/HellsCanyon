#### code/metabolic_analyses/nirK_tree.R ####
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
nirK.tree.unrooted <- read.newick(file = "dataEdited/metabolic_analyses/N/nirK/nirK_for_tree.tree")
nirK.tree.rooted <- midpoint(nirK.tree.unrooted)

#### Dataset with refs: Visualize rooted tree ####
pdf("~/Downloads/nirK_tree.pdf",
    height = 80,
    width = 15)
ggtree(nirK.tree.rooted,
       aes(x = 0, xend = 3)) +
  geom_tiplab(size = 2) +
  geom_text2(aes(subset=!isTip, label=node))
dev.off()



#### Final RAxML tree: Read in tree ####
nirK.raxml.unrooted <- read.newick(file = "dataEdited/metabolic_analyses/N/nirK/RAxML_bipartitions.nirK")
# Final RAxML tree: Rename branches
nirK.raxml.unrooted$tip.label[which(nirK.raxml.unrooted$tip.label %in% names(nirK.metadata.vector))] <- nirK.metadata.vector[nirK.raxml.unrooted$tip.label[which(nirK.raxml.unrooted$tip.label %in% names(nirK.metadata.vector))]]

pdf("~/Downloads/nirK_raxml_raw.pdf",
    height = 100,
    width = 15)
ggtree(nirK.raxml.unrooted,
       aes(x = 0, xend = 3)) +
  geom_tiplab(size = 2) +
  geom_text2(aes(subset=!isTip, label=node))
dev.off()


#### Final RAxML tree: Root tree ####
nirK.tree.rooted <- root(nirK.raxml.unrooted,node = 330)


#### Final RAxML tree: Visualize rooted tree ####
pdf("results/metabolic_analyses/nirK_tree.pdf",
    height = 80,
    width = 15)
ggtree(nirK.tree.rooted,
       aes(x = 0, xend = 3)) +
  geom_tiplab(size = 2) + 
  geom_nodelab(aes(x = branch),
               vjust = -.4,
               hjust = 0.6,
               size = 2)
dev.off()
