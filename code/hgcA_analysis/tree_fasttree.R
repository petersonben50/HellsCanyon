#### code/hgcA_analysis/tree_fasttree.R ####
# Benjamin D. Peterson

# This script will look at our initial FastTree
# of hgcA


#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon/")
library(ggtree)
library(phangorn)
library(readxl)
library(tidyverse)
library(treeio)
cb.translator <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")



#### Check out FastTree of final alignment ####

# Read in tree
tree.name <- "dataEdited/hgcA_analysis/phylogeny/rough_hgcA.tree"
hgcA.tree.unrooted <- read.newick(tree.name)
rm(tree.name)

# Check out unrooted tree
pdf("results/hgcA_analysis/hgcA_tree_FastTree_unrooted.pdf",
    height = 120,
    width = 8)
ggtree(hgcA.tree.unrooted) +
  geom_tiplab(size=2.5, align = TRUE) +
  geom_text2(aes(subset=!isTip, label=node),
             size = 1.5)
dev.off()

# Branch leading to paralogs is 1807



#### Root tree ####
hgcA.tree <- root(hgcA.tree.unrooted,
                  node = 1807,
                  edgelabel = TRUE)



#### List of hgcA sequences ####
hgcA.list <- readLines("dataEdited/hgcA_analysis/hgcA_rep_list.txt")
this.study.indices <- which(hgcA.tree$tip.label %in% hgcA.list)



#### Set color vector ####
color.vector <- rep("black", length(hgcA.tree$tip.label))
color.vector[this.study.indices] <- "red"

# Visualize tree
pdf("results/hgcA_analysis/hgcA_tree_FastTree.pdf",
    height = 120,
    width = 8)
ggtree(hgcA.tree) +
  geom_tiplab(size=2.5,
              align = TRUE,
              col = color.vector) +
  geom_text2(aes(subset = !isTip,
                 label = node),
             size = 1.5)
dev.off()




#### Pull out sequence names to remove ####

nodes.to.remove <- c(1574, 1689, 1555, 1554, 1559, 1402,
                     1382, 1500, 1318, 1344, 1284, 1263,
                     1247, 1355, 1362, 1365, 1176, 1146,
                     1203, 1217, 1241, 1141, 1131, 1102,
                     1086, 1073, 1055, 1702, 1976, 1756,
                     1780, 1723, 1715, 1899, 1922, 1973,
                     1870, 1869, 1819, 1808)
seqs.to.remove <- vector()
for (node.of.interest in nodes.to.remove) {
  tree.subset.to.remove <- tree_subset(hgcA.tree,
                                       node = node.of.interest,
                                       levels_back = 0)
  seqs.to.remove <- c(seqs.to.remove, tree.subset.to.remove$tip.label)
}
fileConn <- file("dataEdited/hgcA_analysis/phylogeny/seqs_to_remove.txt")
writeLines(seqs.to.remove, fileConn)
close(fileConn)

rm(list = ls())
