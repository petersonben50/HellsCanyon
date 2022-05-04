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
tree.name <- "dataEdited/bins/binAnalysis/bacteroidetes_tree/tree_building/RAxML_bipartitions.Bacteroidetes_rp16"
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
hgcA.bin.vector <- readLines("dataEdited/bins/binAnalysis/bacteroidetes_tree/hgcA_bin_list.txt")
naming.vector[which(names(naming.vector) %in% hgcA.bin.vector)] <- paste(naming.vector[which(names(naming.vector) %in% hgcA.bin.vector)],
                                                                         "**",
                                                                         sep = "")



#### Rename tips ####
hgcA.tree$tip.label[hgcA.tree$tip.label %in% names(naming.vector)] <- naming.vector[hgcA.tree$tip.label[hgcA.tree$tip.label %in% names(naming.vector)]]


#### Remove BS values <50 ####
hgcA.tree$node.label <- as.numeric(hgcA.tree$node.label)
hgcA.tree$node.label[hgcA.tree$node.label < 50] <- ""


#### Read in GTDB taxonomy data ####
gtdb.tax.data <- read.table("dataEdited/bins/binAnalysis/bacteroidetes_tree/taxonomy_summary.txt",
                            sep = '\t',
                            header = TRUE)


#### Trim off excess ####
ggtree(hgcA.tree) + 
  geom_tiplab(size=2.5) +
  geom_text2(aes(subset=!isTip, label=node))
# Cut off at node 50
hgcA.tree.subset <- tree_subset(hgcA.tree,
                                node = 50,
                                levels_back = 0)



#### Make color vector for tree ####
color.vector <- rep(cb.translator["bluishgreen"], length(hgcA.tree.subset$tip.label))
color.vector[grep("\\(GCF", hgcA.tree.subset$tip.label)] <- "black"
color.vector[grep("\\(GCA", hgcA.tree.subset$tip.label)] <- "grey50"
color.vector[grep("BAC_00", hgcA.tree.subset$tip.label)] <- cb.translator["skyblue"]
color.vector[grep("anvio", hgcA.tree.subset$tip.label)] <- cb.translator["orange"]



#### Generate tree ####
bin.tree <- ggtree(hgcA.tree.subset, aes(x = 0, xend = 1.75)) + 
  geom_tiplab(size=2,
              align = TRUE,
              colour = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.4,
               hjust = 0.6,
               size = 1) +
  geom_treescale(x = 0.02,
                 y = 21,
                 width = 0.1)
bin.tree


#### Save out tree ###
pdf("results/bins/binAnalysis/phylogeny/bacteroidetes_tree.pdf",
    height = 3.5,
    width = 3.5)
bin.tree
dev.off()
