#### code/binning/phylogenies/bacteroidetes_tree.R ####
# Benjamin D. Peterson

# This script will generate images of
# the phylogenetic trees for the bacteroidetes bins.


#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon/")
library(ggtree)
library(phangorn)
library(readxl)
library(tidyverse)
library(treeio)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")



#### Read in naming info ####
naming.df <- read_xlsx("dataEdited/bins/binAnalysis/phylogeny/bacteroidetes/genome_notes.xlsx")
naming.vector <- paste(naming.df$binName,
                       " (", naming.df$accessionNumber, ")",
                       sep = "")
naming.vector <- gsub(pattern = " \\(NA\\)",
                      replacement = "",
                      naming.vector)
names(naming.vector) <- naming.df$tipLabel





#### FastTree ####

# Read in tree
tree.name <- "dataEdited/bins/binAnalysis/bacteroidetes_tree/tree_building/rp16.tree"
hgcA.tree.unrooted <- read.newick(tree.name)
rm(tree.name)
ggtree(hgcA.tree.unrooted,
       aes(x = 0,
           xend = 1)) + 
  geom_tiplab(size=2.5, align = TRUE)

#### Root tree ####
mrca.flavo <- getMRCA(hgcA.tree.unrooted,
                      c("GCF_000236705.1",
                        "GCF_000194605.1"))
hgcA.tree <- root(hgcA.tree.unrooted,
                  node = mrca.flavo)
rm(mrca.flavo,
   hgcA.tree.unrooted)



#### Rename tips ####
hgcA.tree$tip.label[hgcA.tree$tip.label %in% names(naming.vector)] <- naming.vector[hgcA.tree$tip.label[hgcA.tree$tip.label %in% names(naming.vector)]]



# Check out rooted tree
pdf("dataEdited/bins/binAnalysis/bacteroidetes_tree/bacteroidetes_tree_FastTree_rooted.pdf",
    height = 10,
    width = 5)
ggtree(hgcA.tree,
       aes(x = 0,
           xend = 1)) + 
  geom_tiplab(size=2.5, align = TRUE)
dev.off()
rm(hgcA.tree)



#### Read in RAxML tree ####
# Read in tree
tree.name <- "dataEdited/bins/binAnalysis/phylogeny/bacteroidetes/RAxML_bipartitions.Bacteroidetes_rp16"
hgcA.tree.unrooted <- read.newick(tree.name)
rm(tree.name)


#### Root tree ####
mrca.flavo <- getMRCA(hgcA.tree.unrooted,
                      c("GCF_000236705.1",
                        "GCF_000194605.1"))
hgcA.tree <- root(hgcA.tree.unrooted,
                  node = mrca.flavo)
rm(mrca.flavo)


#### Add hgcA info to renaming vector ####
hgcA.bin.vector <- readLines("dataEdited/bins/binAnalysis/phylogeny/bacteroidetes/hgcA/hgcA_bin_list.txt")
naming.vector[which(names(naming.vector) %in% hgcA.bin.vector)] <- paste(naming.vector[which(names(naming.vector) %in% hgcA.bin.vector)],
                                                                         "**",
                                                                         sep = "")



#### Rename tips ####
hgcA.tree$tip.label[hgcA.tree$tip.label %in% names(naming.vector)] <- naming.vector[hgcA.tree$tip.label[hgcA.tree$tip.label %in% names(naming.vector)]]


#### Make color vector for tree ####
color.vector <- rep(cb.translator["bluishgreen"], length(hgcA.tree$tip.label))
color.vector[grep("\\(GCF", hgcA.tree$tip.label)] <- "black"
color.vector[grep("\\(GCA", hgcA.tree$tip.label)] <- "grey50"
color.vector[grep("BAC_00", hgcA.tree$tip.label)] <- cb.translator["skyblue"]
color.vector[grep("anvio", hgcA.tree$tip.label)] <- cb.translator["orange"]


#### Remove BS values <50 ####
hgcA.tree$node.label <- as.numeric(hgcA.tree$node.label)
hgcA.tree$node.label[hgcA.tree$node.label < 50] <- ""



#### Generate tree ####
bin.tree <- ggtree(hgcA.tree, aes(x = 0, xend = 1.25)) + 
  geom_tiplab(size=2.5,
              align = TRUE,
              colour = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.4,
               hjust = 0.6,
               size = 2) +
  geom_treescale(x = 0.05,
                 y = 30,
                 width = 0.2)

#### Save out tree ###
pdf("results/bins/binAnalysis/phylogeny/bacteroidetes_tree.pdf",
    height = 6,
    width = 5)
bin.tree
dev.off()
