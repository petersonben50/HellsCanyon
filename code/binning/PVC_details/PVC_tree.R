#### code/binning/PVC_details/PVC_tree.R ####
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



#### Investigate GTDB tree ####
GTDB.tree <- read.newick("dataEdited/bins/binAnalysis/phylogeny/gtdbtk.bac120.classify.tree")
mrca.node <- getMRCA(phy = GTDB.tree,
                     c("anvio_hgcA_0261", "anvio_hgcA_0220"))

# Subset tree
PVC.gtdb.tree <- tree_subset(GTDB.tree,
                             node = mrca.node,
                             levels_back = 0)
ggtree(PVC.gtdb.tree,
       aes(x = 0,
           xend = 2))  + 
  geom_tiplab(size=2.5, align = TRUE)
















#### Initial check of tree ####

# Read in tree
tree.name <- "dataEdited/binning/phylogeny/PVC/RAxML_bipartitions.PVC_rp16"
hgcA.tree.unrooted <- read.newick(tree.name)
rm(tree.name)

# Check out unrooted tree
pdf("dataEdited/binning/phylogeny/PVC/PVC_tree_FastTree_unrooted.pdf",
    height = 10,
    width = 5)
ggtree(hgcA.tree.unrooted,
       aes(x = 0,
           xend = 2)) + 
  geom_tiplab(size=2.5, align = TRUE)
dev.off()



#### Read in naming info ####
naming.df <- read.table("dataEdited/binning/phylogeny/PVC/tip_naming.tsv",
                        sep = '\t',
                        header = TRUE)
naming.vector <- paste(naming.df$binName,
                       " (", naming.df$accessionNumber, ")",
                       sep = "")
naming.vector <- gsub(pattern = " \\(NA\\)",
                      replacement = "",
                      naming.vector)
names(naming.vector) <- naming.df$tipLabel



#### Root tree ####
mrca.plancto <- getMRCA(hgcA.tree.unrooted,
                        c("GCF_000255705.1",
                          "GCF_000025185.1"))
hgcA.tree <- root(hgcA.tree.unrooted,
                  node = mrca.plancto)


#### Find Lentispharae node ####
lent.node <- getMRCA(hgcA.tree,
                     naming.df %>%
                       filter(binPhylum == "Lentisphaerae") %>%
                       select(tipLabel) %>%
                       unlist(use.names = FALSE))
kir.node <- getMRCA(hgcA.tree,
                     naming.df %>%
                       filter(binPhylum == "Kiritimatiellaeota") %>%
                       select(tipLabel) %>%
                       unlist(use.names = FALSE))


# Rename tips
hgcA.tree$tip.label[hgcA.tree$tip.label %in% names(naming.vector)] <- naming.vector[hgcA.tree$tip.label[hgcA.tree$tip.label %in% names(naming.vector)]]


#### Make color vector for tree ####
color.vector <- rep(cb.translator["bluishgreen"], length(hgcA.tree$tip.label))
color.vector[grep("\\(GCF", hgcA.tree$tip.label)] <- "black"
color.vector[grep("\\(GCA", hgcA.tree$tip.label)] <- "grey50"
color.vector[grep("KIR|LEN", hgcA.tree$tip.label)] <- cb.translator["skyblue"]
color.vector[grep("anvio", hgcA.tree$tip.label)] <- cb.translator["orange"]


#### Remove BS values <50 ####
hgcA.tree$node.label <- as.numeric(hgcA.tree$node.label)
hgcA.tree$node.label[hgcA.tree$node.label < 50] <- ""






#### Generate tree ####
bin.tree <- ggtree(hgcA.tree) + 
  geom_tiplab(size=2.5,
              align = TRUE,
              col = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 2) +
  geom_treescale(x = 0.05,
                 y = 30,
                 width = 0.25) +
  geom_cladelabel(node = lent.node, label="Lentisphaerae",
                  align = T, geom = 'label') +
  geom_cladelabel(node = kir.node, label="Kiritimatiellaeota",
                  align = T, geom = 'label')

pdf("results/binning/phylogeny/PVC_tree_unedited.pdf",
    height = 10,
    width = 5)
bin.tree
dev.off()
