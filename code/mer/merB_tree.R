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


#### Check out FastTree of alignment of original sequences ####
tree.name <- "dataEdited/mer/identification/merB_tree_for_ID.tree"
merB.tree.unrooted <- read.newick(tree.name)
rm(tree.name)


#### Root tree ####
merB.tree <- midpoint(merB.tree.unrooted, node.labels = "support")


#### Color vector ####
color.vector <- rep("black", length(merB.tree$tip.label))
color.vector[c(grep("fall2017", merB.tree$tip.label),
               grep("HC18", merB.tree$tip.label),
               grep("KMBP00", merB.tree$tip.label))] <- "red"



#### Visualize tree ####
pdf("results/mer/merB_tree_identification.pdf",
    height = 75,
    width = 8)
merB.tree.image <- ggtree(merB.tree) + 
  geom_tiplab(size=2.5,
              align = TRUE,
              col = color.vector) +
  geom_text2(aes(subset = !isTip,
                 label = node),
             size = 1.5)
merB.tree.image
dev.off()


#### Pull out sequence lists ####

# Divergent cluster 1: node = 1307
DC1_tree <- tree_subset(merB.tree,
                        node = 1307,
                        levels_back = 0)
writeLines(DC1_tree$tip.label,
           "dataEdited/mer/identification/DC1_list.txt")

# merB cluster 1: node = 1421
merB1_tree <- tree_subset(merB.tree,
                        node = 1421,
                        levels_back = 0)
writeLines(merB1_tree$tip.label,
           "dataEdited/mer/identification/merB1_list.txt")


# merB cluster 2: node = 1331
merB2_tree <- tree_subset(merB.tree,
                          node = 1331,
                          levels_back = 0)
ggtree(merB2_tree,
       aes(x = 0,
           xend = 25)) + 
  geom_tiplab(size=2.5,
              align = TRUE)
writeLines(merB2_tree$tip.label,
           "dataEdited/mer/identification/merB2_list.txt")


# merB cluster 3: node = 1029
merB3_tree <- tree_subset(merB.tree,
                          node = 1029,
                          levels_back = 0)
writeLines(merB3_tree$tip.label,
           "dataEdited/mer/identification/merB3_list.txt")





#### Check out FastTree of alignment of original sequences ####
tree.name <- "dataEdited/mer/identification/subclusters/all_seqs_of_interest.tree"
merB.subset.tree.unrooted <- read.newick(tree.name)
rm(tree.name)

#### Root tree ####
merB.subset.tree <- midpoint(merB.subset.tree.unrooted, node.labels = "support")


#### Color vector ####
color.vector <- rep("black", length(merB.subset.tree$tip.label))
color.vector[c(grep("fall2017", merB.subset.tree$tip.label),
               grep("HC18", merB.subset.tree$tip.label),
               grep("KMBP00", merB.subset.tree$tip.label))] <- "red"

#### Generate a tree image ####
merB.tree.image <- ggtree(merB.subset.tree) + 
  geom_tiplab(size=2.5,
              align = TRUE,
              col = color.vector) +
  geom_text2(aes(subset = !isTip,
                 label = node),
             size = 1.5)

pdf("results/mer/merB_tree_subset_with_alignment.pdf",
    height = 30,
    width = 8)
msaplot(merB.tree.image,
        "dataEdited/mer/identification/subclusters/all_seqs_of_interest_trimmed.afa",
        offset = 8)
dev.off()
