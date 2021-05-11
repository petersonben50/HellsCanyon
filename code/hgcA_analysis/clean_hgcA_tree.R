#### code/hgcA_analysis/clean_hgcA_tree.R ####
# Benjamin D. Peterson

# This script will look at our phylogenetic tree 
# iterations and ultimately clean up the final one.


#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon/")
library(ggtree)
library(phangorn)
library(readxl)
library(tidyverse)
library(treeio)
cb.translator <- readRDS("~/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")



#### Check out FastTree of final alignment ####

# Read in tree
tree.name <- "dataEdited/hgcA_analysis/phylogeny/rough_hgcA.tree"
hgcA.tree.unrooted <- read.newick(tree.name)
rm(tree.name)

# Check out unrooted tree
pdf("results/hgcA_analysis/hgcA_tree_FastTree_unrooted.pdf",
    height = 120,
    width = 8)
ggtree(hgcA.tree.unrooted) + 
  geom_tiplab(size=2.5, align = TRUE) +
  geom_text2(aes(subset=!isTip, label=node),
             size = 1.5)
dev.off()

# Branch leading to paralogs is 1807



#### Root tree ####
hgcA.tree <- root(hgcA.tree.unrooted,
                  node = 1807,
                  edgelabel = TRUE)



#### List of hgcA sequences ####
hgcA.list <- readLines("dataEdited/hgcA_analysis/hgcA_rep_list.txt")
this.study.indices <- which(hgcA.tree$tip.label %in% hgcA.list)



#### Set color vector ####
color.vector <- rep("black", length(hgcA.tree$tip.label))
color.vector[this.study.indices] <- "red"

# Visualize tree
pdf("results/hgcA_analysis/hgcA_tree_FastTree.pdf",
    height = 120,
    width = 8)
ggtree(hgcA.tree) + 
  geom_tiplab(size=2.5,
              align = TRUE,
              col = color.vector) +
  geom_text2(aes(subset = !isTip,
                 label = node),
             size = 1.5)
dev.off()




#### Pull out sequence names to remove ####

nodes.to.remove <- c(1574, 1689, 1555, 1554, 1559, 1402,
                     1382, 1500, 1318, 1344, 1284, 1263,
                     1247, 1355, 1362, 1365, 1176, 1146,
                     1203, 1217, 1241, 1141, 1131, 1102,
                     1086, 1073, 1055, 1702, 1976, 1756,
                     1780, 1723, 1715, 1899, 1922, 1973,
                     1870, 1869, 1819, 1808)
seqs.to.remove <- vector()
for (node.of.interest in nodes.to.remove) {
  tree.subset.to.remove <- tree_subset(hgcA.tree,
                                       node = node.of.interest,
                                       levels_back = 0)
  seqs.to.remove <- c(seqs.to.remove, tree.subset.to.remove$tip.label)
}
fileConn <- file("dataEdited/hgcA_analysis/phylogeny/seqs_to_remove.txt")
writeLines(seqs.to.remove, fileConn)
close(fileConn)

rm(list = ls())






##################################################
#### Clean up RAxML tree from final alignment ####
##################################################

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
tree.name <- "dataEdited/hgcA_analysis/phylogeny/RAxML_bipartitions.hgcA"
hgcA.tree.unrooted <- read.newick(tree.name)
rm(tree.name)


#### Fix tip labels ####
hgcA.tree.unrooted$tip.label <- paste(hgcA.tree.unrooted$tip.label %>% strsplit("_") %>% sapply("[", 1),
                                      hgcA.tree.unrooted$tip.label %>% strsplit("_") %>% sapply("[", 2),
                                      hgcA.tree.unrooted$tip.label %>% strsplit("_") %>% sapply("[", 3),
                                      sep = "_")

# Root tree
hgcA.tree <- root(phy = hgcA.tree.unrooted,
                  outgroup = c("paralog_Thermosulfurimonas_dismutans",
                               "paralog_Candidatus_Omnitrophica"),
                  edgelabel = TRUE)


# Get indices
reference.indices.mate <- which(hgcA.tree$tip.label %in% names(MATE.renaming.vector))
this.study.indices <- which(hgcA.tree$tip.label %in% hgcA.list)

# Change names:
hgcA.tree$tip.label[reference.indices.mate] <- MATE.renaming.vector[hgcA.tree$tip.label[reference.indices.mate]]
hgcA.tree$tip.label[this.study.indices] <- hgcA.renaming.vector[hgcA.tree$tip.label[this.study.indices]]


# Set color vector
color.vector <- rep("grey", length(hgcA.tree$tip.label))
color.vector[reference.indices.mate] <- colorblind.color.vector["black"]
color.vector[this.study.indices] <- colorblind.color.vector["vermillion"]

# Visualize tree
pdf("results/hgcA_analysis/hgcA_tree_RAxML_rooted.pdf",
    height = 60,
    width = 10)
ggtree(hgcA.tree, aes(x = 0,
                      xend = 9)) + 
  geom_tiplab(size=2.5, align = TRUE, col = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 2)
dev.off()




#### Save out tree ####
saveRDS(hgcA.tree,
        "dataEdited/hgcA_analysis/phylogeny/hgcA_clean_tree.rds")


#### Save out color vector ####
saveRDS(color.vector,
        "dataEdited/hgcA_analysis/phylogeny/hgcA_clean_tree_color_vector.rds")
