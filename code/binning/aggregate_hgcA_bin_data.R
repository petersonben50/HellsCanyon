#### code/binning/aggregate_hgcA_bin_data.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(igraph)
library(tidyverse)
library(vegan)



#### Calculate HMSs ####

ANI.values <- read.table(file = "dataEdited/bins/binning/bins_hgcA/binsHgcA.all.ani.out.cleaned",
                         sep = "\t",
                         header = TRUE,
                         stringsAsFactors = FALSE) %>%
  mutate(GENOME1 = gsub(".fna", "", GENOME1),
         GENOME2 = gsub(".fna", "", GENOME2))
names(ANI.values)[3:6] <- c("ANI_1", "ANI_2", "AF_1", "AF_2")

# Check plot of values
plot(y = ANI.values$ANI_1,
     x = ANI.values$AF_1,
     xlim = c(0.2, 1.0),
     ylim = c(70, 100),
     pch = 18)
points(y = ANI.values$ANI_2,
       x = ANI.values$AF_2,
       pch = 18)
abline(h = 97)
abline(v = 0.5)


# Set cut-offs
ANI.cutoff <- 97
cov.cutoff <- 0.5

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



#### Add in assembly and year information to each bin ####
autoBins <- grep(pattern = "hgcA",
                 x = allData$binID,
                 invert = TRUE,
                 value = TRUE)
temp.binAssemblyInfo <- data.frame(strsplit(autoBins, "_") %>% sapply("[", 1),
                                   autoBins)
names(temp.binAssemblyInfo) <- c("assemblyID", "binID")
binAssemblyInfo <- rbind(read.table("dataEdited/bins/binning/manualBinning/notes/bin_to_assembly.tsv",
                                    col.names = c("assemblyID", "binID")),
                         temp.binAssemblyInfo)
rm(autoBins, temp.binAssemblyInfo)
yearData <- read.csv("metadata/lists/assembly_key.csv",
                     col.names = c("assemblyID", "year"))

allData <- left_join(allData,
                     binAssemblyInfo) %>%
  left_join(yearData)
rm(binAssemblyInfo,
   yearData)


#### Read in taxonomy data ####
taxonomy.data <- read.table("dataEdited/bins/binning/bins_hgcA/taxonomy_summary.txt",
                            sep = "\t",
                            header = FALSE,
                            col.names = c("binID", "gtdb_tax"),
                            stringsAsFactors = FALSE)
allData <- allData %>%
  left_join(taxonomy.data)
rm(taxonomy.data)


#### CheckM data ####

checkm.data <- read.csv("dataEdited/bins/binning/bins_hgcA/checkM_stats.csv",
                        stringsAsFactors = FALSE)
names(checkm.data) <- c("binID",
                        paste("checkM_",
                              c("completeness", "contamination", "strain_het",
                                "genome_length", "scaf_count", "N50", "mean_scaf_length",
                                "longest_scaf", "GC", "predicted_genes"),
                              sep = ""))
allData <- allData %>%
  left_join(checkm.data)
rm(checkm.data)



#### Save out data ####
write.csv(allData,
          "dataEdited/bins/binning/bins_hgcA/bin_dereplication_data.csv",
          row.names = FALSE,
          quote = FALSE)
