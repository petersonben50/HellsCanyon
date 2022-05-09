#### code/binning/PVC_details/hgcA_tree.R ####
# Benjamin D. Peterson

# This script will generate images of
# the phylogenetic trees for the PVC bins.


#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon/")
library(ggtree)
library(phangorn)
library(readxl)
library(tidyverse)
library(treeio)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")



#### Read in hgcA tree ####
hgcA.tree.file <- "dataEdited/bins/binAnalysis/PVC_details/hgcA_analysis/tree_generation/RAxML_bipartitions.PVC_hgcA"
hgcA.tree.unrooted <- read.newick(hgcA.tree.file)
hgcA.tree <- midpoint(hgcA.tree.unrooted)



#### Read in G2B file ####
g2b.file <- read.table("dataEdited/bins/binAnalysis/PVC_details/hgcA_analysis/PVC_hgcA_G2B.tsv",
                       sep = '\t',
                       col.names = c("geneID", "binID"))
g2b.vector <- g2b.file$binID
names(g2b.vector) <- g2b.file$geneID



#### Read in naming info ####
naming.df <- read_xlsx("dataEdited/bins/binAnalysis/PVC_details/tree_generation/tip_naming.xlsx") %>%
  filter(tipLabel %in% g2b.vector)
naming.vector <- paste(naming.df$binName,
                       " (", naming.df$accessionNumber, ")",
                       sep = "")
naming.vector <- gsub(pattern = " \\(NA\\)",
                      replacement = "",
                      naming.vector)
names(naming.vector) <- naming.df$tipLabel



#### Rename tips ####
tips.to.be.renamed <- which(hgcA.tree$tip.label %in% names(g2b.vector))
hgcA.tree$tip.label[tips.to.be.renamed] <- naming.vector[g2b.vector[hgcA.tree$tip.label[tips.to.be.renamed]]]



#### Make color vector for tree ####
color.vector <- rep(cb.translator["bluishgreen"], length(hgcA.tree$tip.label))
color.vector[grep("\\*\\*", hgcA.tree$tip.label)] <- cb.translator["skyblue"]
color.vector[grep("\\(GCF", hgcA.tree$tip.label)] <- "black"
color.vector[grep("\\(GCA", hgcA.tree$tip.label)] <- "grey50"
color.vector[grep("KIR|LEN", hgcA.tree$tip.label)] <- cb.translator["orange"]
color.vector[grep("IMG:", hgcA.tree$tip.label)] <- cb.translator["vermillion"]
names(color.vector) <- NULL


#### Remove BS values <50 ####
hgcA.tree$node.label <- as.numeric(hgcA.tree$node.label)
hgcA.tree$node.label[hgcA.tree$node.label < 50] <- ""

#### Generate tree ####
hgcA.tree.visualized <- ggtree(hgcA.tree,
                               aes(x = 0,
                                   xend = 2.5)) + 
  geom_tiplab(size=2.5,
              align = TRUE,
              colour = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 1.8) +
  geom_treescale(x = 0.05,
                 y = 20,
                 width = 0.25)
hgcA.tree.visualized



#### Save out tree ####
pdf("results/bins/binAnalysis/PVC_details/hgcA_tree.pdf",
    height = 4,
    width = 3.5)
hgcA.tree.visualized
dev.off()
