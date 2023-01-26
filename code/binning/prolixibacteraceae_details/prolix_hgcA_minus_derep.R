#### code/binning/ANI/prolix_derep.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(igraph)
library(readxl)
library(tidyverse)
library(vegan)


#### Read in taxonomy data ####
prolix.tax.data <- read_xlsx("dataEdited/bins/binning/autoBinning/taxonomy_summary.xlsx") %>%
  mutate(family = classification %>%
           strsplit(";f__") %>% sapply("[", 2) %>%
           strsplit(";g__") %>% sapply("[", 1),
         genus = classification %>%
           strsplit(";g__") %>% sapply("[", 2) %>%
           strsplit(";s__") %>% sapply("[", 1)) %>%
  filter(family == "Prolixibacteraceae") %>%
  rename(binID = user_genome)


#### Read in taxonomy data ####
prolix.quality.data <- read.csv("dataEdited/bins/binning/autoBinning/checkM_stats.csv") %>%
  rename(binID = Bin.Id,
         completeness = Completeness,
         redundancy = Contamination,
         genomeSize = Genome.size..bp.,
         geneCount = X..predicted.genes) %>%
  select(binID, completeness, redundancy, genomeSize, geneCount)


#### Calculate HMSs ####
ANI.values <- read.table(file = "dataEdited/bins/binning/autoBinning/hqBins.all.ani.out.cleaned",
                         sep = "\t",
                         header = TRUE,
                         stringsAsFactors = FALSE) %>%
  mutate(GENOME1 = gsub(".fna", "", GENOME1),
         GENOME2 = gsub(".fna", "", GENOME2)) %>%
  filter(GENOME1 %in% prolix.tax.data$binID &
           GENOME2 %in% prolix.tax.data$binID)
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


#### Combine data ####
allData <- allData %>%
  left_join(prolix.quality.data) %>%
  left_join(prolix.tax.data)




#### Bins to keep ####
bins.to.keep.phylogeny <- c("HC18HY02_bin_0055", "KMBP004F_bin_0124", "KMBP004F_bin_0265",
                            "fall2017coassembly_bin_0080", "KMBP004F_bin_0183", "fall2017cluster6_bin_0136",
                            "fall2017coassembly_bin_0181", "fall2017coassembly_bin_0389")
# Leaving out HC18ME02_bin_0076 because it's the same as the anvio_hgcA_0130
writeLines(bins.to.keep.phylogeny,
           "dataEdited/bins/binAnalysis/phylogeny/prolixibacteraceae/hgcA_minus_bins_to_use_for_phylogeny.txt")



#### Dereplicate depth data ####
depth.data <- readRDS("dataEdited/bins/binAnalysis/hqBins/bin_depth_clean.rds") %>%
  filter(binID %in% allData$binID)

allData.with.depth <- allData %>%
  select(HMS, binID, classification) %>%
  right_join(depth.data) %>%
  group_by(HMS, metagenomeID) %>%
  summarise(coverage = mean(coverage))
# Add binID of bins used in phylogeny
key.df <- allData %>%
  filter(binID %in% bins.to.keep.phylogeny) %>%
  select(binID, HMS)
allData.with.depth <- allData.with.depth %>%
  left_join(key.df)

#### Save out data ####
saveRDS(allData.with.depth,
        "dataEdited/bins/binAnalysis/phylogeny/prolixibacteraceae/depth_of_bins.rds")
