#### code/binning/taxonomy.R ####
# Benjamin D. Peterson


#### Clean up crew on line 5 ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(ggtree)
library(phangorn)
library(tidyverse)
library(treeio)


#### Read in tree ####
GTDB.tree <- read.newick("dataEdited/binning/taxonomy/gtdbtk.bac120.classify.tree")
binNames <- grep("anvio_hgcA",
                 GTDB.tree$tip.label,
                 value = TRUE)


#### Isolate GTDB tree for Bacteroidetes bin ####
bin <- "anvio_hgcA_0130"
tree_file_name <- paste("dataEdited/binning/taxonomy/",
                        bin,
                        "_gtdbtk_tree.pdf",
                        sep = "")
subset <- GTDB.tree %>%
  tree_subset(bin,
              levels_back = 4) %>%
  midpoint(node.labels = "support")
color.vector.to.use <- rep("black",
                           length(subset$tip.label))
color.vector.to.use[grep(bin, subset$tip.label)] <- "red"

subset %>%
  ggtree() +
  geom_tiplab(size = 2,
              colour = color.vector.to.use) +
  geom_text2(aes(subset=!isTip, label=node))

subset <- subset %>%
  tree_subset(node = 1079,
              levels_back = 0)
color.vector.to.use <- rep("black", 
                           length(subset$tip.label))
color.vector.to.use[grep(bin, subset$tip.label)] <- "red"

print(subset %>%
  ggtree() +
  geom_tiplab(size = 2,
              colour = color.vector.to.use) +
  geom_text2(aes(subset=!isTip, label=node)))
write(subset$tip.label %>%
        gsub("RS_", "", .) %>%
        gsub("GB_", "", .),
      file = "dataEdited/binning/taxonomy/bacteroidetes_references.txt")
