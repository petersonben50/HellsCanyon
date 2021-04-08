#### code/binning/aggregate_bin_data.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(igraph)
library(tidyverse)
library(vegan)



#### Calculate HMSs ####

ANI.values <- read.table(file = "dataEdited/binning/ANI/goodBins.all.ani.out.cleaned",
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
  filter(ANI_1 > ANI.cutoff &
           ANI_2 > ANI.cutoff &
           AF_1 > cov.cutoff &
           AF_2 > cov.cutoff) %>%
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
HMS.bin.info <- do.call("rbind",
                        list(HMS.ID, lone.bins))
rm(ANI.cutoff, cov.cutoff, HMS.ID, lone.bins,
   edgelist, ANI.values, adjacency.graph,
   bin.list)




#### Add in assembly and year information to each bin ####

binAssemblyInfo <- read.table("dataEdited/binning/manualBinning/notes/bin_to_assembly.tsv",
                              col.names = c("assemblyID", "binID"))
yearData <- read.csv("metadata/lists/assembly_key.csv",
                     col.names = c("assemblyID", "year"))
allData <- left_join(HMS.bin.info,
                     binAssemblyInfo) %>%
  left_join(yearData)
rm(binAssemblyInfo,
   yearData,
   HMS.bin.info)


#### Read in taxonomy data ####

taxonomy.data <- read.table("dataEdited/binning/taxonomy/gtdbtk.bac120.summary.tsv",
                            sep = "\t",
                            header = TRUE,
                            stringsAsFactors = FALSE) %>%
  rename(binID = user_genome,
         gtdb_tax = classification) %>%
  select(binID, gtdb_tax)
allData <- allData %>%
  left_join(taxonomy.data)
rm(taxonomy.data)


#### CheckM data ####

checkm.data <- read.csv("dataEdited/binning/quality/checkM_stats.csv",
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


#### anvi'o completeness/redundancy ####

anvio.quality <- read.table("dataEdited/binning/quality/bins_summary_hgcA.txt",
                            sep = "\t",
                            stringsAsFactors = FALSE,
                            header = TRUE)
names(anvio.quality) <- paste("anvio_",
                              names(anvio.quality),
                              sep = "")
anvio.quality <- anvio.quality %>%
  rename(binID = anvio_bins)
allData <- allData %>%
  left_join(anvio.quality)
rm(anvio.quality)


#### Save out data ####

write.csv(allData,
          "dataEdited/binning/bin_data.csv",
          row.names = FALSE,
          quote = FALSE)
