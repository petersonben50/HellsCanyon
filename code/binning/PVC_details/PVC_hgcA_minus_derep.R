#### code/binning/PVC_details/PVC_hgcA_minus_derep.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(igraph)
library(readxl)
library(tidyverse)
library(vegan)


#### Calculate HMSs ####
ANI.values <- read.table(file = "dataEdited/bins/binAnalysis/PVC_details/genome_dereplication/bins_to_dereplicate.all.ani.out.cleaned",
                         sep = "\t",
                         header = TRUE,
                         stringsAsFactors = FALSE) %>%
  mutate(GENOME1 = gsub(".fna", "", GENOME1),
         GENOME2 = gsub(".fna", "", GENOME2))
names(ANI.values)[3:6] <- c("ANI_1", "ANI_2", "AF_1", "AF_2")

# Check plot of values
plot(y = ANI.values$ANI_1,
     x = ANI.values$AF_1,
     xlim = c(0.1, 1.0),
     ylim = c(70, 100),
     pch = 18)
points(y = ANI.values$ANI_2,
       x = ANI.values$AF_2,
       pch = 18)
abline(h = 97)
abline(v = 0.4)




# Set cut-offs
ANI.cutoff <- 97
cov.cutoff <- 0.4

# Generate edgelist
edgelist <- ANI.values %>%
  filter((ANI_1 > ANI.cutoff |
            ANI_2 > ANI.cutoff) &
           (AF_1 > cov.cutoff |
              AF_2 > cov.cutoff)) %>%
  select(GENOME1, GENOME2)
# Generate the graph from the edgelist
adjacency.graph <- graph_from_data_frame(edgelist)
# Store the bin and HMS name in a df.
HMS.ID <- data.frame(paste("HMS.",
                           clusters(adjacency.graph)$membership,
                           sep = ""),
                     names(clusters(adjacency.graph)$membership),
                     stringsAsFactors = FALSE)
names(HMS.ID) <- c("HMS", "binID")
HMS.ID <- HMS.ID %>%
  mutate(binID = binID %>%
           strsplit(".fna") %>%
           sapply("[", 1)) %>%
  arrange(HMS)


# Check out loner bins
bin.list <- unique(c(ANI.values$GENOME1, ANI.values$GENOME2)) %>%
  strsplit(".fna") %>%
  sapply("[", 1)
lone.bins <- data.frame(bin.list[!(bin.list %in% HMS.ID$bin)],
                        bin.list[!(bin.list %in% HMS.ID$bin)],
                        stringsAsFactors = FALSE)
names(lone.bins) <- c("HMS", "binID")
# Add loner bins to list
allData <- do.call("rbind",
                   list(HMS.ID, lone.bins))
rm(ANI.cutoff, cov.cutoff, HMS.ID, lone.bins,
   edgelist, ANI.values, adjacency.graph,
   bin.list)


#### Keep hgcA+ bins ####
bins.to.keep.phylogeny <- allData %>%
  filter(grepl("anvio", binID)) %>%
  select(binID) %>%
  unlist(use.names = FALSE)


#### Remove HMSs with an hgcA+ bin from remaining analysis ####
hgcA.HMS <- allData %>%
  filter(grepl("anvio", binID)) %>%
  select(HMS) %>%
  unlist(use.names = FALSE)
allData <- allData %>%
  filter(!(HMS %in% hgcA.HMS))



#### Read in taxonomy data ####
PVC.tax.data <- read_xlsx("dataEdited/bins/binning/autoBinning/taxonomy_summary.xlsx") %>%
  mutate(phylum = classification %>%
           strsplit(";p__") %>% sapply("[", 2) %>%
           strsplit(";c__") %>% sapply("[", 1),
         class = classification %>%
           strsplit(";c__") %>% sapply("[", 2) %>%
           strsplit(";o__") %>% sapply("[", 1)) %>%
  filter(class %in% c("Kiritimatiellae", "Lentisphaeria")) %>%
  rename(binID = user_genome)


#### Read in checkM data ####
PVC.quality.data <- read.csv("dataEdited/bins/binning/autoBinning/checkM_stats.csv") %>%
  rename(binID = Bin.Id,
         completeness = Completeness,
         redundancy = Contamination,
         genomeSize = Genome.size..bp.,
         geneCount = X..predicted.genes) %>%
  select(binID, completeness, redundancy, genomeSize, geneCount)


#### Combine data ####
allData <- allData %>%
  left_join(PVC.quality.data) %>%
  left_join(PVC.tax.data)




#### Bins to keep ####
bins.to.keep.phylogeny <- c(bins.to.keep.phylogeny, "HC18HY300_bin_0036", "HC18HY300_bin_0043", "HC18HY300_bin_0050",
                            "fall2017cluster6_bin_0092", "fall2017cluster6_bin_0034", "fall2017cluster6_bin_0083",
                            "HC18HY300_bin_0028", "HC18HY300_bin_0070", "HC18ME02_bin_0081",
                            "HC18ME02_bin_0105", "fall2017coassembly_bin_0332")
writeLines(bins.to.keep.phylogeny,
           "dataEdited/bins/binAnalysis/PVC_details/genome_dereplication/bins_to_keep.txt")



# #### Dereplicate depth data ####
# depth.data <- readRDS("dataEdited/bins/binAnalysis/hqBins/bin_depth_clean.rds") %>%
#   filter(binID %in% allData$binID)
# 
# allData.with.depth <- allData %>%
#   select(HMS, binID, classification) %>%
#   right_join(depth.data) %>%
#   group_by(HMS, metagenomeID) %>%
#   summarise(coverage = mean(coverage))
# saveRDS(allData.with.depth,
#         "dataEdited/bins/binAnalysis/phylogeny/PVCibacteraceae/depth_of_bins.rds")
