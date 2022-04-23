#### code/binning/prolixibacteraceae_details/prolixibacteraceae_hgcA_tree.R ####
# Benjamin D. Peterson

# This script will generate images of
# the phylogenetic trees for the hgcA
# sequences from the prolixibacteraceae bins.


#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon/")
library(ggtree)
library(phangorn)
library(readxl)
library(tidyverse)
library(treeio)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")



#### Prepare vector for tip labels ####
gene.to.bin <- read.table("dataEdited/bins/binAnalysis/prolixibacteraceae_details/hgcA/hgcA_G2B.tsv",
                          sep = '\t',
                          col.names = c("geneID", "binID"))
tax.data.bact <- read.table("dataEdited/bins/binAnalysis/phylogeny/bacteroidetes/taxonomy_summary.txt",
                            sep = '\t',
                            header = TRUE) %>%
  filter(user_genome %in% gene.to.bin$binID) %>%
  mutate(family = classification %>%
           strsplit("f__") %>% sapply("[", 2) %>%
           strsplit(";g__") %>% sapply("[", 1)) %>%
  rename(binID = user_genome) %>%
  mutate(accessionID = binID) %>%
  select(binID, family, accessionID)

# Tax data for prolixibacteraceae
tax.data.prol <- read_xlsx("dataEdited/bins/binAnalysis/phylogeny/prolixibacteraceae/taxonomy_summary.xlsx") %>%
  filter(user_genome %in% gene.to.bin$binID) %>%
  rename(binID = user_genome) %>%
  mutate(accessionID = binID) %>%
  select(binID, family, accessionID)

# Taxonomic data for Mendota study sequences
tax.data.5M <- read.table("dataEdited/bins/binAnalysis/prolixibacteraceae_details/hgcA/gtdb_taxonomy_summary.tsv",
                          sep = '\t',
                          header = TRUE) %>%
  filter(binID %in% grep("BAC_", binID, value = TRUE)) %>%
  select(binID, family) %>%
  mutate(family = family %>% gsub("f__", "", .)) %>%
  mutate(accessionID = "Mendota_study")

# Taxonomic data for Mendota study sequences
tax.data.this.study <- read.table("dataEdited/bins/binning/bins_hgcA/taxonomy_summary.txt",
                          sep = '\t',
                          col.names = c("binID", "classification")) %>%
  mutate(family = classification %>%
           strsplit("f__") %>% sapply("[", 2) %>%
           strsplit(";g__") %>% sapply("[", 1)) %>%
  filter(binID == "anvio_hgcA_0130") %>%
  select(binID, family) %>%
  mutate(accessionID = "This_study")

all.tax.data <- rbind(tax.data.this.study, tax.data.5M, tax.data.prol, tax.data.bact)
rm(tax.data.this.study, tax.data.5M, tax.data.prol, tax.data.bact)


all.tax.data <- all.tax.data %>%
  left_join(gene.to.bin)

# Prepare renaming vector
renaming.vector <- paste(all.tax.data$family, " (", all.tax.data$binID, ")",
                         sep = '')
names(renaming.vector) <- all.tax.data$geneID



#### Investigate FastTree tree ####
hgcA.bact.tree.unrooted <- read.newick("dataEdited/bins/binAnalysis/prolixibacteraceae_details/hgcA/hgcA_bact.tree")
# Check out unrooted tree
ggtree(hgcA.bact.tree.unrooted,
       aes(x = 0,
           xend = 0.4)) +
  geom_tiplab(size = 2)


#### Root tree ####
fastTree.tree <- midpoint(hgcA.bact.tree.unrooted)


#### Rename branches ####
fastTree.tree$tip.label <- renaming.vector[fastTree.tree$tip.label]


#### Visualize Prolixibacteraceae tree ####
ggtree(fastTree.tree,
       aes(x = 0,
           xend = 1.2)) +
  geom_tiplab(size = 2)





#### Read in RAxML tree ####
raxmlTree.tree.unrooted <- read.newick("dataEdited/bins/binAnalysis/prolixibacteraceae_details/hgcA/RAxML_bipartitions.hgcA_bact")


#### Root tree ####
raxmlTree.tree <- midpoint(raxmlTree.tree.unrooted)



#### Rename branches ####
raxmlTree.tree$tip.label <- renaming.vector[raxmlTree.tree$tip.label]



#### Set up color vector ####
color.vector <- rep(cb.translator["bluishgreen"], length(raxmlTree.tree$tip.label))
color.vector[grep("GCF_", raxmlTree.tree$tip.label)] <- "black"
color.vector[grep("GCA_", raxmlTree.tree$tip.label)] <- "grey50"
color.vector[grep("anvio", raxmlTree.tree$tip.label)] <- cb.translator["skyblue"]


#### Save out tree ####
pdf("results/bins/binAnalysis/prolixibacteraceae_details/hgcA_tree.pdf",
    height = 11,
    width = 8.5)
ggtree(raxmlTree.tree,
       aes(x = 0,
           xend = 1.8)) +
  geom_tiplab(size = 2,
              colour = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.4,
               hjust = 0.6,
               size = 2) +
  geom_treescale(x = 0.0,
                 y = 22,
                 width = 0.2)
dev.off()
