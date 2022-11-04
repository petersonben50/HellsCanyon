#### code/binning/geobacterales/geobacterales_tree.R ####
# Benjamin D. Peterson

# This script will generate images of
# the phylogenetic trees for the Geobacterales bins.


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
geo.tree <- tree_subset(GTDB.tree,
                        "anvio_hgcA_0210",
                        levels_back = 6)
# Save out tree
pdf("dataEdited/bins/binAnalysis/phylogeny/geobacterales/GTDB_tree.pdf",
    width = 6,
    height = 6)
ggtree(geo.tree,
       aes(x = 0,
           xend = 0.45)) +
  geom_tiplab(size = 2)
dev.off()
geo.tree$tip.label


#### Pull out accession numbers of reference bins ####

writeLines(c(tree_subset(GTDB.tree,
                         "anvio_hgcA_0210",
                         levels_back = 5)$tip.label %>%
               gsub("GB_", "", .) %>%
               gsub("RS_", "", .),
             "GCF_000016745.1",
             "GCF_000007985.2"),
           "dataEdited/bins/binAnalysis/phylogeny/geobacterales/GTDB_ref_accession_numbers.txt")




#### Read in naming info ####
naming.df <- read_xlsx("dataEdited/bins/binAnalysis/phylogeny/geobacterales/genome_notes_geo_ref.xlsx")
naming.vector <- paste(naming.df$species,
                       " (", naming.df$user_genome, ")",
                       sep = "")
names(naming.vector) <- naming.df$user_genome


#### Add hgcA info to renaming vector ####
hgcA.bin.vector <- readLines("dataEdited/bins/binAnalysis/phylogeny/geobacterales/hgcA/hgcA_bin_list.txt")
naming.vector[which(names(naming.vector) %in% hgcA.bin.vector)] <- paste(naming.vector[which(names(naming.vector) %in% hgcA.bin.vector)],
                                                                         "**",
                                                                         sep = "")


#### Check out FastTree ####
FastTree.unrooted <- read.newick("dataEdited/bins/binAnalysis/phylogeny/geobacterales/rp16.tree")
mrca.flavo <- getMRCA(FastTree.unrooted,
                      c("GCF_000016745.1",
                        "GCF_000007985.2"))
FastTree <- root(FastTree.unrooted,
                 node = mrca.flavo)

#### Rename tip labels ####
ref.vector.index <- which(FastTree$tip.label %in% names(naming.vector))
FastTree$tip.label[ref.vector.index] <- naming.vector[FastTree$tip.label[ref.vector.index]]

#### Visualize tree ####
ggtree(FastTree,
       aes(x = 0, xend = 0.6)) + 
  geom_tiplab(size=2.5,
              align = TRUE)




#### Check out RAxML tree ####
raxml.unrooted <- read.newick("dataEdited/bins/binAnalysis/phylogeny/geobacterales/RAxML_bipartitions.geobacterales_rp16")
mrca.flavo <- getMRCA(raxml.unrooted,
                      c("GCF_000016745.1",
                        "GCF_000007985.2"))
raxml.tree <- root(raxml.unrooted,
                   node = mrca.flavo)

#### Rename tip labels ####
ref.vector.index <- which(raxml.tree$tip.label %in% names(naming.vector))
raxml.tree$tip.label[ref.vector.index] <- naming.vector[raxml.tree$tip.label[ref.vector.index]]


#### Remove BS values <50 ####
raxml.tree$node.label <- as.numeric(raxml.tree$node.label)
raxml.tree$node.label[raxml.tree$node.label < 50] <- ""


#### Visualize and save out tree ####
pdf("results/bins/binAnalysis/geobacterales/geobacterales_tree.pdf",
    height = 4,
    width = 5)
ggtree(raxml.tree,
       aes(x = 0, xend = 1.5)) + 
  geom_tiplab(size=2.5,
              align = TRUE) + 
  geom_nodelab(aes(x = branch),
               vjust = -.4,
               hjust = 0.6,
               size = 2.5) +
  geom_treescale(x = 0.0,
                 y = 15,
                 width = 0.1)

dev.off()
