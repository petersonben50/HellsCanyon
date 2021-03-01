#### code/mash/NMDS_plot_MGs_2017.R ####
# Benjamin D. Peterson

#### Always start with a clean slate ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(tidyverse)
library(vegan)




#### Read in distance data ####
dist.data.2017 <- read.table("dataEdited/metagenomes/mash_data/HCC_MG_2017.dist",
                             stringsAsFactors = FALSE,
                             sep = "\t")
names(dist.data.2017) <- c("refID", "queryID", "mash_dist",
                           "pvalue", "matching_hashes")
clean.dist.data.2017 <- dist.data.2017 %>%
  mutate(refID = refID %>%
           gsub("temp_MG_files/", "", .) %>%
           gsub(".fastq.gz", "", .)) %>%
  mutate(queryID = queryID %>%
           gsub("temp_MG_files/", "", .) %>%
           gsub(".fastq.gz", "", .)) %>%
  select(refID, queryID, mash_dist) %>%
  spread(key = queryID,
         value = mash_dist)



#### Process metadata ####
MG.data <- read.csv("metadata/metagenome_metadata.csv",
                    stringsAsFactors = FALSE) %>%
  mutate(RM = gsub("RM", "", RM))
MG.data.vector <- paste("RM", MG.data$RM, ",",
                        MG.data$depth, "m",
                        sep = "")
names(MG.data.vector) <- MG.data$metagenomeID





#### Generate ordination ####

# Prepare matrix
clean.dist.data.2017.matrix <- clean.dist.data.2017 %>%
  select(-refID)
row.names(clean.dist.data.2017.matrix) <- clean.dist.data.2017$refID

# Run NMDS function
nmds.sites <- metaMDS(clean.dist.data.2017.matrix)

# Set up data to plot
data.scores = as.data.frame(scores(nmds.sites),
                            stringsAsFactors = FALSE)
data.scores$site.info <- paste(rownames(data.scores), '\n',
                               MG.data.vector[rownames(data.scores)],
                               sep = "")
data.scores$RM <- MG.data.vector[rownames(data.scores)] %>%
  strsplit(",") %>%
  sapply("[", 1)


#### Make plot ####

pdf("results/metagenomes/mash/2017_NMDS_MGs.pdf",
    height = 6,
    width = 6)
ggplot(data.scores, aes(x = NMDS1, y = NMDS2, col = RM)) + 
  geom_point(size = 0) + 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + 
  labs(x = "NMDS1", y = "NMDS2") +
  geom_text(aes(label = site.info)) +
  ggtitle("NMDS based on Jaccard distance\nof 2017 metagenomes (by Mash)")
dev.off()
