#### code/binning/metabolism/MoORs.R ####
# Benjamin D. Peterson


#### Get set up ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon/")
# library(ape)
library(dplyr)
library(ggtree)
library(phangorn)
library(readxl)
library(treeio)
source("/Users/benjaminpeterson/Documents/programs/R/functions/reading_HMM_output.R")


#### Generate bin/gene conversion vectors ####
G2B <- read.table("dataEdited/bins/binAnalysis/bin_ORFs/ORFs_G2B.tsv",
                  stringsAsFactors = FALSE,
                  sep = "\t",
                  header = FALSE)
names(G2B) <- c("geneID", "binID")
# Gene to bin vector
gene.to.bin.vector <- G2B$binID
names(gene.to.bin.vector) <- G2B$geneID

# Bin to gene vector
bin.to.gene.vector <- G2B$geneID
names(bin.to.gene.vector) <- G2B$binID

rm(G2B)


#### Read in tree ####
MoOR.tree <- read.newick("dataEdited/bins/binAnalysis/metabolism/MoORs/putative_MoORs_1.tree")


#### Read in list of bins to use ####
binList <- readLines("dataEdited/bins/binning/bins_hgcA_keepers/bins_hgcA_keepers_list.txt")


#### Read in HMM data ####
HMM.data <- read_tblout("dataEdited/bins/binAnalysis/metabolism/MoORs/MoOR.out") %>%
  select(domainName, seq_score)
hmm.score.vector <- HMM.data$seq_score
names(hmm.score.vector) <- HMM.data$domainName
rm(HMM.data)



#### Read in taxonomy data ####
bin.data <- read_xlsx("dataEdited/bins/binning/bins_hgcA/bin_dereplication_data_edits.xlsx",
                      sheet = "bin_dereplication_data_trimmed") %>%
  select(binID, HMS, gtdb_tax) %>%
  filter(binID %in% binList)

bin.taxonomy.vector <- strsplit(bin.data$gtdb_tax,
                                split = ";") %>% sapply("[", 3)
names(bin.taxonomy.vector) <- bin.data$binID
rm(bin.data)


#### Read in HMM outputs ####
HMM.output <- read_tblout(file = "dataEdited/bins/binAnalysis/metabolism/MoORs/MoOR.out") %>%
  select(domain_name, sequence_score) %>%
  rename(seqID = domain_name,
         HMM_score = sequence_score) %>%
  as.data.frame()



#### Generate metadata vector ####
metadata <- read.table("references/MoORs/MoOR_key.tsv",
                       stringsAsFactors = FALSE,
                       sep = '\t',
                       header = FALSE)
names(metadata) <- c("geneID", "name")

ref.ID.to.name.vector <- make.unique(metadata$name)
names(ref.ID.to.name.vector) <- metadata$geneID

ref.name.to.ID.vector <- make.unique(metadata$geneID)
names(ref.name.to.ID.vector) <- metadata$name



#### Generate renaming vector ####
renaming.vector <- MoOR.tree$tip.label
names(renaming.vector) <- MoOR.tree$tip.label

# Renaming for my bins
my.bin.index <- which(renaming.vector %in% bin.to.gene.vector)
renaming.vector[my.bin.index] <- paste(renaming.vector[my.bin.index],
                                       " (", hmm.score.vector[renaming.vector[my.bin.index]], ")",
                                       " - ",
                                       gene.to.bin.vector[renaming.vector[my.bin.index]],
                                       " (", bin.taxonomy.vector[gene.to.bin.vector[renaming.vector[my.bin.index]]], ")",
                                       sep = "")

# Renaming for references
reference.index <- which(renaming.vector %in% names(ref.ID.to.name.vector))
renaming.vector[reference.index] <- ref.ID.to.name.vector[renaming.vector[reference.index]]
rm(reference.index)



#### Make color vector ####
color.vector <- rep("black", length(MoOR.tree$tip.label))
color.vector[my.bin.index] <- "red"



#### Write out new tree with updated labels ####
MoOR.tree.for.printing <- MoOR.tree
MoOR.tree.for.printing$tip.label <- renaming.vector[MoOR.tree$tip.label]
write.tree(MoOR.tree.for.printing,
           "dataEdited/bins/binAnalysis/metabolism/MoORs/FastTree_MoOR_edited.tree")



#### Inspect this tree ####
pdf("dataEdited/bins/binAnalysis/metabolism/MoORs/MoOR_tree_1.pdf",
    height = 30,
    width = 8)
ggtree(MoOR.tree.for.printing) + 
  geom_tiplab(size = 2,
              col = color.vector) +
  geom_text2(aes(subset=!isTip,
                 label = node)) +
  xlim(c(0, 10))
dev.off()
# Node 207 leads to large cluster of hits with no close references
# These hits also have low scores.
# We'll remove these from further analysis.


#### Identify sequences to remove ####
to.remove.tree <- tree_subset(MoOR.tree.for.printing,
                              207,
                              levels_back = 0)
to.remove.list <- names(to.remove.tree$tip.label)
MoORs.list <-  readLines("dataEdited/bins/binAnalysis/metabolism/MoORs/putative_MoORs_list.txt")
to.keep.list <- MoORs.list[!(MoORs.list %in% to.remove.list)]
writeLines(text = to.keep.list,
           con = "dataEdited/bins/binAnalysis/metabolism/MoORs/MoORs_list.txt") 



#### Read in RAxML tree ####
MoOR.tree.unrooted <- read.newick("dataEdited/bins/binAnalysis/metabolism/MoORs/tree_take_2/RAxML_bipartitions.MoORs")
MoOR.tree <- midpoint(MoOR.tree.unrooted)
MoOR.tree.for.printing <- MoOR.tree
MoOR.tree.for.printing$tip.label <- renaming.vector[MoOR.tree$tip.label]


#### Make bin index ####
my.bin.index <- which(renaming.vector %in% bin.to.gene.vector)



#### Make color vector ####
color.vector <- rep("black", length(MoOR.tree$tip.label))
color.vector[my.bin.index] <- "red"


# Check it out in Geneious 
pdf("dataEdited/bins/binAnalysis/metabolism/MoORs/tree_take_2/MoOR_RAxML_tree.pdf",
    height = 30,
    width = 8)
ggtree(MoOR.tree.for.printing,
       aes(x = 0, xend = 10)) + 
  geom_tiplab(size = 2,
              col = color.vector) + 
  geom_nodelab(aes(x = branch),
               vjust = -.3,
               size = 2)
dev.off()

