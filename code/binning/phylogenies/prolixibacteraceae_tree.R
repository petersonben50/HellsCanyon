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


#### Set up color vector ####
color.vector <- rep(cb.translator["bluishgreen"], length(fastTree.tree$tip.label))
color.vector[grep("GCF_", fastTree.tree$tip.label)] <- "black"
color.vector[grep("GCA_", fastTree.tree$tip.label)] <- "grey50"
color.vector[grep("\\*\\*", fastTree.tree$tip.label)] <- cb.translator["reddishpurple"]
color.vector[grep("anvio", fastTree.tree$tip.label)] <- cb.translator["skyblue"]



#### Visualize Prolixibacteraceae tree ####
ggtree(fastTree.tree,
       aes(x = 0,
           xend = 0.4)) +
  geom_tiplab(size = 2,
              colour = color.vector)





#### Read in RAxML tree ####
raxmlTree.tree.unrooted <- read.newick("dataEdited/bins/binAnalysis/phylogeny/prolixibacteraceae/RAxML_bipartitions.prolixibacteraceae_rp16")


#### Root tree ####
mrca.flavo <- getMRCA(raxmlTree.tree.unrooted,
                      c("GCF_000236705.1",
                        "GCF_000194605.1"))
raxmlTree.tree <- root(raxmlTree.tree.unrooted,
                      node = mrca.flavo)


#### Read in taxonomy info ####
tax.info <- read_xlsx("dataEdited/bins/binAnalysis/phylogeny/prolixibacteraceae/taxonomy_summary.xlsx")
tax.info.vector <- paste(tax.info$family, "-", tax.info$genus, " (", tax.info$user_genome, ")",
                         sep = "")
names(tax.info.vector) <- tax.info$user_genome


#### Add hgcA information ####
hgcA.data <- readLines("dataEdited/bins/binAnalysis/phylogeny/prolixibacteraceae/hgcA/hgcA_bin_list.txt")
hgcA.index <- which(names(tax.info.vector) %in% hgcA.data)
tax.info.vector[hgcA.index] <- paste(tax.info.vector[hgcA.index], "**", sep = "")


# Rename references
ref.index <- which(raxmlTree.tree$tip.label %in% tax.info$user_genome)
raxmlTree.tree$tip.label[ref.index] <- tax.info.vector[raxmlTree.tree$tip.label[ref.index]]



#### Set up color vector ####
color.vector <- rep(cb.translator["bluishgreen"], length(raxmlTree.tree$tip.label))
color.vector[grep("GCF_", raxmlTree.tree$tip.label)] <- "black"
color.vector[grep("GCA_", raxmlTree.tree$tip.label)] <- "grey50"
color.vector[grep("\\*\\*", raxmlTree.tree$tip.label)] <- cb.translator["reddishpurple"]
color.vector[grep("anvio", raxmlTree.tree$tip.label)] <- cb.translator["skyblue"]


#### Save out tree ####
pdf("results/bins/binAnalysis/phylogeny/prolixibacteraceae_tree.pdf",
    height = 11,
    width = 8.5)
ggtree(raxmlTree.tree,
       aes(x = 0,
           xend = 1.2)) +
  geom_tiplab(size = 2,
              colour = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.4,
               hjust = 0.6,
               size = 2) +
  geom_treescale(x = 0.05,
                 y = 100,
                 width = 0.2)
dev.off()
