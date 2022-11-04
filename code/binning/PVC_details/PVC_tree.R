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



####-------------Investigate GTDB tree-------------####
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

rm(GTDB.tree, PVC.gtdb.tree, mrca.node)






####-------------Visualize ML tree from RAxML-------------####


#### Read in PVC bin tree ####
PVC.bin.tree.unrooted <- read.newick("dataEdited/bins/binAnalysis/PVC_details/tree_generation/RAxML_bipartitions.HCC_PVC")


#### Read in list of hgcA+ bins ####
hgcA.G2B.table <- read.table("dataEdited/bins/binAnalysis/PVC_details/hgcA_analysis/PVC_hgcA_G2B.tsv",
                             sep = '\t',
                             col.names = c("geneID", "binID"))



#### Read in naming info ####
naming.df <- read_xlsx("dataEdited/bins/binAnalysis/PVC_details/tree_generation/tip_naming.xlsx") %>%
  filter(tipLabel %in% PVC.bin.tree.unrooted$tip.label)
naming.vector <- paste(naming.df$binName,
                       " (", naming.df$accessionNumber, ")",
                       sep = "")
naming.vector <- gsub(pattern = " \\(NA\\)",
                      replacement = "",
                      naming.vector)
names(naming.vector) <- naming.df$tipLabel

index.of.hgcA.in.naming.vector <- which(names(naming.vector) %in% hgcA.G2B.table$binID)
naming.vector[index.of.hgcA.in.naming.vector] <- paste(naming.vector[index.of.hgcA.in.naming.vector], "**",
                                                       sep = "")



#### Root tree ####
mrca.plancto <- getMRCA(PVC.bin.tree.unrooted,
                        c("GCF_000255705.1",
                          "GCF_000025185.1"))
PVC.bin.tree <- root(PVC.bin.tree.unrooted,
                     node = mrca.plancto)


#### Find Lentispharae node ####
lent.node <- getMRCA(PVC.bin.tree,
                     naming.df %>%
                       filter(binPhylum == "Lentisphaerae") %>%
                       select(tipLabel) %>%
                       unlist(use.names = FALSE))
kir.node <- getMRCA(PVC.bin.tree,
                     naming.df %>%
                       filter(binPhylum == "Kiritimatiellaeota") %>%
                       select(tipLabel) %>%
                       unlist(use.names = FALSE))


# Rename tips
PVC.bin.tree$tip.label[PVC.bin.tree$tip.label %in% names(naming.vector)] <- naming.vector[PVC.bin.tree$tip.label[PVC.bin.tree$tip.label %in% names(naming.vector)]]


#### Make color vector for tree ####
color.vector <- rep(cb.translator["bluishgreen"], length(PVC.bin.tree$tip.label))
color.vector[grep("\\*\\*", PVC.bin.tree$tip.label)] <- cb.translator["vermillion"]
color.vector[grep("\\(GCF", PVC.bin.tree$tip.label)] <- "black"
color.vector[grep("\\(GCA", PVC.bin.tree$tip.label)] <- "grey50"
color.vector[grep("KIR|LEN", PVC.bin.tree$tip.label)] <- cb.translator["skyblue"]
color.vector[grep("IMG:", PVC.bin.tree$tip.label)] <- cb.translator["blue"]


#### Remove BS values <50 ####
high.bootstrap.index <- which(as.numeric(PVC.bin.tree$node.label) >= 80)
low.bootstrap.index <- which(as.numeric(PVC.bin.tree$node.label) < 80)
PVC.bin.tree$node.label[high.bootstrap.index] <- "*"
PVC.bin.tree$node.label[low.bootstrap.index] <- ""



#### Generate tree ####
bin.tree <- ggtree(PVC.bin.tree,
                   aes(x = 0,
                       xend = 4)) + 
  geom_tiplab(size=2.5,
              align = TRUE,
              colour = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = 0.2,
               size = 3.5) +
  geom_treescale(x = 0.05,
                 y = 30,
                 width = 0.25) +
  geom_cladelabel(node = lent.node, label="Lentisphaerae",
                  align = T,
                  offset = 3.6,
                  offset.text	= 0.15,
                  hjust = 0.5,
                  angle = 270) +
  geom_cladelabel(node = kir.node, label="Kiritimatiellaeota",
                  align = T,
                  offset = 3.6,
                  offset.text	= 0.15,
                  hjust = 0.5,
                  angle = 270)
bin.tree



#### Save tree to PDF ####
pdf("results/bins/binAnalysis/PVC_details/PVC_tree.pdf",
    height = 6,
    width = 3.5)
bin.tree
dev.off()
