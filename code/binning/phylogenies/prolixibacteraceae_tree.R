#### code/binning/phylogenies/prolixibacter_tree.R ####
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


#### Investigate GTDB tree ####
GTDB.tree <- read.newick("dataEdited/bins/binAnalysis/phylogeny/gtdbtk.bac120.classify.tree")

# Subset tree
prolix.tree <- tree_subset(GTDB.tree,
                         "anvio_hgcA_0130",
                         levels_back = 2)
ggtree(prolix.tree,
       aes(x = 0,
           xend = 0.4)) +
  geom_tiplab(size = 2)
prolix.tree$tip.label


#### Read in FastTree ####
fastTree.tree.unrooted <- read.newick("dataEdited/bins/binAnalysis/phylogeny/prolixibacteraceae/rp16.tree")


#### Root tree ####
mrca.flavo <- getMRCA(fastTree.tree.unrooted,
                      c("GCF_000236705.1",
                        "GCF_000194605.1"))
fastTree.tree <- root(fastTree.tree.unrooted,
                      node = mrca.flavo)
ggtree(fastTree.tree,
       aes(x = 0,
           xend = 0.4)) +
  geom_tiplab(size = 2,
              colour = color.vector)


#### Read in RAxML tree ####
fastTree.tree.unrooted <- read.newick("dataEdited/bins/binAnalysis/phylogeny/prolixibacteraceae/RAxML_bipartitions.prolixibacteraceae_rp16")


#### Root tree ####
mrca.flavo <- getMRCA(fastTree.tree.unrooted,
                      c("GCF_000236705.1",
                        "GCF_000194605.1"))
fastTree.tree <- root(fastTree.tree.unrooted,
                      node = mrca.flavo)


#### Read in taxonomy info ####
tax.info <- read_xlsx("dataEdited/bins/binAnalysis/phylogeny/prolixibacteraceae/taxonomy_summary.xlsx")
tax.info.vector <- paste(tax.info$family, tax.info$genus, sep = "-")
names(tax.info.vector) <- tax.info$user_genome


#### Add hgcA information ####
hgcA.data <- readLines("dataEdited/bins/binAnalysis/phylogeny/prolixibacteraceae/hgcA/hgcA_bin_list.txt")
hgcA.index <- which(names(tax.info.vector) %in% hgcA.data)
tax.info.vector[hgcA.index] <- paste(tax.info.vector[hgcA.index], "**", sep = "")


# Rename references
ref.index <- which(fastTree.tree$tip.label %in% tax.info$user_genome)
fastTree.tree$tip.label[ref.index] <- tax.info.vector[fastTree.tree$tip.label[ref.index]]


#### Set up color vector ####
color.vector <- rep("black", length(fastTree.tree$tip.label))
color.vector[grep("\\*\\*", fastTree.tree$tip.label)] <- "red"
color.vector[grep("anvio", fastTree.tree$tip.label)] <- "blue"

ggtree(fastTree.tree,
       aes(x = 0,
           xend = 0.4)) +
  geom_tiplab(size = 2,
              colour = color.vector)
