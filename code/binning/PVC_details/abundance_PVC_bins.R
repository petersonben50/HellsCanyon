#### code/binning/PVC_details/abundance_PVC_bins.R ####
# Benjamin D. Peterson

# This script will plot the abundance of all the
# PVC bins from this study.

#### Always start with a clean slate ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon/")
library(lubridate)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")
cb.translator <- c(cb.translator, "grey50")
names(cb.translator)[length(cb.translator)] <- "grey"



####-------------Read in MHC data-------------####
high.mhc.bins <- readRDS("dataEdited/bins/binAnalysis/PVC_details/metabolism/high_MHC_bins.rds")
low.mhc.bins <- readRDS("dataEdited/bins/binAnalysis/PVC_details/metabolism/low_MHC_bins.rds")
MHC.content <- c(rep("high", length(high.mhc.bins)),
                 rep("low", length(low.mhc.bins)))
names(MHC.content) <- c(high.mhc.bins, low.mhc.bins)
rm(high.mhc.bins, low.mhc.bins)


####-------------Read in taxonomy depth data-------------####
tax.data <- read.table("dataEdited/bins/binAnalysis/PVC_details/taxonomy/taxonomy_summary.txt",
                       sep = '\t', header = TRUE) %>%
  rename(binID = user_genome) %>%
  mutate(class = classification %>%
           strsplit("__") %>% sapply("[", 4) %>%
           strsplit(";") %>% sapply("[", 1),
         family = classification %>%
           strsplit("__") %>% sapply("[", 5) %>%
           strsplit(";") %>% sapply("[", 1)) %>%
  select(binID, class, family)


####-------------Read in hgcA data-------------####
hgcA.list <- readLines("dataEdited/bins/binning/bins_hgcA_keepers/bins_hgcA_keepers_list.txt")




####-------------Read in depth data-------------####
depth.data <- rbind(readRDS("dataEdited/bins/binAnalysis/hqBins/bin_depth_clean.rds"),
                    readRDS("dataEdited/bins/binAnalysis/depth/bin_depth_clean.rds")) %>%
  filter(binID %in% names(MHC.content)) %>%
  mutate(mhcContent = MHC.content[binID])



####-------------Add in metadata-------------####
MG.metadata <- read_xlsx("metadata/metagenome_metadata.xlsx")
all.data <- depth.data %>%
  left_join(MG.metadata) %>%
  left_join(tax.data) %>%
  mutate(hgcA = (binID %in% hgcA.list))



####-------------Set up color vector-------------####
color.vector <- cb.translator[c("vermillion", "yellow", "bluishgreen")]
names(color.vector) <- c("Kiritimatiellae-high", "Kiritimatiellae-low", "Lentisphaeria-low")
shape.vector <- c(15, 4)
names(shape.vector) <- c(TRUE, FALSE)


####-------------Plot out data-------------####
all.data %>%
  filter(((RM %in% c(286, 300)) & (year(date) %in% c(2017, 2018))) |
           ((RM %in% c(300, 310)) & (year(date) == 2019)))%>%
  ggplot(aes(y = coverage,
             x = depth,
             group = binID)) +
  geom_point(aes(color = paste(class, mhcContent, sep = "-"),
                 shape = hgcA)) +
  geom_line(aes(color = paste(class, mhcContent, sep = "-"))) +
  scale_color_manual(values = color.vector,
                     name = "Phylum - MHC content") +
  scale_shape_manual(values = shape.vector,
                     name = "hgcA present") +
  facet_wrap(~year(date) + RM, nrow = 3) +
  scale_y_continuous(trans = "log10",
                     limits = c(0.0001, 3.5)) +
  coord_flip(xlim = c(80, 0)) +
  theme_bw()

