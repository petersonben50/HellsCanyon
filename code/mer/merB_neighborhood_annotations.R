#### code/mer/merB_neighborhood_annotations.R ####
# Written by Benjamin D. Peterson


#### Get ready. Get set ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(lubridate)
library(readxl)
library(tidyverse)


#### Read in annotations ####
merB.neighborhood.annotations <- read.table("dataEdited/mer/scaffolds/friendly_neighborhood_genes_annotationsAll.tsv",
                                            skip = 1,
                                            sep = '\t')
names(merB.neighborhood.annotations) <- c("dump", "seqID", "KOFAM", "threshold", "score", "evalue", "geneInfo")
merB.neighborhood.annotations.filtered <- merB.neighborhood.annotations %>%
  group_by(seqID) %>%
  slice_max(order_by = score, n = 1) %>%
  mutate(annotationInfo = paste(geneInfo, "-",
                                KOFAM,
                                sep = "")) %>%
  select(seqID, annotationInfo)


#### Read in GFF files ####
merB.gff <- read.table("dataEdited/mer/scaffolds/merB_geneNeighborhood_raw.gff")
names(merB.gff) <- c("scaffoldID", "source", "type", "start", "end",
                     "score", "strandedness", "zero", "annotation")
merB.gff <- merB.gff %>%
  mutate(seqID = paste(scaffoldID, "_",
                       annotation %>%
                         strsplit(";") %>% sapply("[", 1) %>%
                         strsplit("_") %>% sapply("[", 2),
                       sep = "")) %>%
  left_join(merB.neighborhood.annotations.filtered) %>%
  mutate(annotation = paste(annotation, "product=",
                            annotationInfo,
                            sep = "")) %>%
  select(scaffoldID, source, type, start, end,
         score, strandedness, zero, annotation)
merB.gff$annotation[grep("product=NA$", merB.gff$annotation)] <- paste(merB.gff$annotation[grep("product=NA$", merB.gff$annotation)],
                                                                      "-NA",
                                                                      sep = "")

#### Write out GFF table ####
write.table(merB.gff,
            "dataEdited/mer/scaffolds/merB_geneNeighborhood_annot.gff",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
