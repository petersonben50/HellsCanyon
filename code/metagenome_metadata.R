#### code/metagenome_metadata.R ####
# Benjamin D. Peterson

# This file contains the code to generate nice tables with
# information on the metagenome sequencing and assembly.

#### Always start with a clean slate, ya bum ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(dplyr)
library(gridExtra)
library(readxl)


#### Generate mapping key ####
MG.list <- read.csv("metadata/lists/metagenome_key.csv",
                    header = FALSE,
                    col.names = c("metagenomeID", "samplingYear"))
assembly.list <- read.csv("metadata/lists/assembly_key.csv",
                          header = FALSE,
                          col.names = c("assemblyID", "samplingYear"))
mapping.key <- full_join(MG.list,
                         assembly.list) %>%
  select(metagenomeID, assemblyID, samplingYear)

write.table(mapping.key,
            "metadata/lists/mapping_key.tsv",
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)


#### Generate table of MG sample site information ####

# Read in metagenome metadata
MG.metadata <- read_xlsx("dataEdited/dnaSequencing/metagenome_key.xlsx")

# Read in extraction metadata
extraction.metadata <- read_xlsx("dataEdited/dnaExtractions/DNA_extractions_data.xlsx")

# Read in filter metadata
filter.metadata <- read_xlsx("metadata/NA_filters.xlsx")

all.metadata <- left_join(MG.metadata,
                          extraction.metadata,
                          by = "extractionID") %>%
  left_join(filter.metadata,
            by = "filterID")



#### Save out metadata for analysis ####

all.metadata %>%
  select(metagenomeID, date, RM, depth) %>%
  mutate(RM = gsub("RM", "", RM))%>%
  write.csv("metadata/metagenome_metadata.csv",
            row.names = FALSE)



