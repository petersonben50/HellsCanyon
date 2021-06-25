#### code/hgcA_analysis/clean_hgcA_tree.R ####
# Benjamin D. Peterson

# This script will visualize our final hgcA tree


#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon/")
library(ggtree)
library(phangorn)
library(readxl)
library(tidyverse)
library(treeio)
cb.translator <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Read in needed files ####
hgcA.list <- readLines("dataEdited/hgcA_analysis/hgcA_rep_list.txt")
colorblind.color.vector <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Read in Hg-MATE seq info ####
hg.mate.metadata <- read_xlsx("/Users/benjaminpeterson/Documents/research/Hg_MATE/versions/v1.01142021/Hg-MATE-Db.v1.01142021_catalogue.xlsx") %>%
  rename(group = `microbial group`,
         name = `[Organism Name]_Phylum-Class`) %>%
  mutate(name = gsub("sp._", "sp.", name)) %>%
  mutate(name = paste(name %>% strsplit("_") %>% sapply("[", 1),
                      name %>% strsplit("_") %>% sapply("[", 2),
                      sep = "_"),
         treeID = paste(group,
                        " (", name, "-", MATE_id, ")",
                        sep = ""))
MATE.renaming.vector <- hg.mate.metadata$treeID
names(MATE.renaming.vector) <- hg.mate.metadata$MATE_id



#### Read in info on seqs from Jones et al, 2019 ####
jones.metadata <- read.table("references/jones_bin_names.tsv",
                             sep = '\t',
                             col.names = c("binID", "binName")) %>%
  full_join(read.table("references/jones_genes2bin.tsv",
                       sep = '\t',
                       col.names = c("geneID", "binID"))) %>%
  mutate(binName = paste(gsub("Unclassified", "", binName),
                         sep = ""))
jones.renaming.vector <- jones.metadata$binName
names(jones.renaming.vector) <- jones.metadata$geneID


#### Read in sequence data ####
hgcA.df <- read_xlsx("dataEdited/hgcA_analysis/hgcA_dereplication.xlsx")
hgcA.renaming.vector <- paste(hgcA.df$seqID, " (cluster: ", hgcA.df$clstr, ")",
                              sep = "")
names(hgcA.renaming.vector) <- hgcA.df$seqID


#### Read in hgcA to bin info ####
hgcA.to.bin <- read.table("dataEdited/binning/hgcA/hgcA_to_bin.tsv",
                          sep = '\t',
                          header = TRUE)
hgcA.to.bin.vector <- hgcA.to.bin$binID
names(hgcA.to.bin.vector) <- hgcA.to.bin$hgcA_ID
binned.hgcA.index <- which(names(hgcA.renaming.vector) %in% names(hgcA.to.bin.vector))
hgcA.renaming.vector[binned.hgcA.index] <- paste(hgcA.renaming.vector[binned.hgcA.index], " - ",
                                                 hgcA.to.bin.vector[names(hgcA.renaming.vector)[binned.hgcA.index]])


#### Read in tree ####
tree.name <- "dataEdited/hgcA_analysis/phylogeny/final/RAxML_bipartitions.hgcA"
hgcA.tree.unrooted <- read.newick(tree.name)
rm(tree.name)


#### Fix tip labels ####
hgcA.tree.unrooted$tip.label <- paste(hgcA.tree.unrooted$tip.label %>% strsplit("_") %>% sapply("[", 1),
                                      hgcA.tree.unrooted$tip.label %>% strsplit("_") %>% sapply("[", 2),
                                      hgcA.tree.unrooted$tip.label %>% strsplit("_") %>% sapply("[", 3),
                                      sep = "_")
hgcA.tree.unrooted$tip.label <- gsub("_NA_NA", "", hgcA.tree.unrooted$tip.label)

# Root tree
hgcA.tree <- root(phy = hgcA.tree.unrooted,
                  outgroup = c("paralog_Thermosulfurimonas_dismutans",
                               "paralog_Candidatus_Omnitrophica"),
                  edgelabel = TRUE)


# Get indices
reference.indices.mate <- which(hgcA.tree$tip.label %in% names(MATE.renaming.vector))
this.study.indices <- which(hgcA.tree$tip.label %in% hgcA.list)
reference.indices.jones <- which(hgcA.tree$tip.label %in% names(jones.renaming.vector))

# Set color vector
color.vector <- rep(colorblind.color.vector["reddishpurple"], length(hgcA.tree$tip.label))
color.vector[reference.indices.mate] <- colorblind.color.vector["black"]
color.vector[this.study.indices] <- colorblind.color.vector["vermillion"]
color.vector[reference.indices.jones] <- cb.translator["skyblue"]

# Visualize tree
pdf("results/hgcA_analysis/hgcA_tree_RAxML_rooted.pdf",
    height = 18,
    width = 10)
ggtree(hgcA.tree, aes(x = 0,
                      xend = 9)) +
  geom_tiplab(size=2.5, align = TRUE,
              colour = color.vector) +
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 2)
dev.off()

# Change names:
hgcA.tree$tip.label[reference.indices.mate] <- MATE.renaming.vector[hgcA.tree$tip.label[reference.indices.mate]]
hgcA.tree$tip.label[this.study.indices] <- hgcA.renaming.vector[hgcA.tree$tip.label[this.study.indices]]
hgcA.tree$tip.label[reference.indices.jones] <- jones.renaming.vector[hgcA.tree$tip.label[reference.indices.jones]]


#### Save out tree ####
saveRDS(hgcA.tree,
        "dataEdited/hgcA_analysis/phylogeny/hgcA_clean_tree.rds")


#### Save out color vector ####
saveRDS(color.vector,
        "dataEdited/hgcA_analysis/phylogeny/hgcA_clean_tree_color_vector.rds")



#### Remove root ####
hgcA.tree.rootless <- drop.tip(hgcA.tree,
                               grep("paralog", hgcA.tree$tip.label))


reference.indices.mate <- which(hgcA.tree.rootless$tip.label %in% MATE.renaming.vector)
reference.indices.jones <- which(hgcA.tree.rootless$tip.label %in% jones.renaming.vector)
this.study.indices <- which(hgcA.tree.rootless$tip.label %in% hgcA.renaming.vector)

#### Set color vector ####
color.vector.to.use <- rep(cb.translator["reddishpurple"], length(hgcA.tree.rootless$tip.label))
color.vector.to.use[reference.indices.mate] <- cb.translator["vermillion"]
color.vector.to.use[this.study.indices] <- cb.translator["black"]
color.vector.to.use[reference.indices.jones] <- cb.translator["skyblue"]




#### Make image of tree without root ####
pdf("results/hgcA_analysis/hgcA_tree_RAxML_rootRemoved.pdf",
    height = 11,
    width = 8.5)
ggtree(hgcA.tree.rootless,
       aes(x = 0,
           xend = 7)) +
  geom_tiplab(size = 1.9,
              align = TRUE,
              colour = color.vector.to.use) +
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 1.8) +
  geom_treescale(x = 0.05,
                 y = 80,
                 width = 0.5)
dev.off()

