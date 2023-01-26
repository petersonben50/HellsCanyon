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
  unique()


####-------------Read in geochem data to get elevation data-------------####
elevation.info <- read.csv("dataEdited/geochem/geochem_WC.csv") %>%
  select(date, RM, depth, elevation_m) %>%
  mutate(RM = as.character(RM)) %>%
  unique()


####-------------Add in metadata-------------####
MG.metadata <- read_xlsx("metadata/metagenome_metadata.xlsx")
all.data <- depth.data %>%
  left_join(MG.metadata) %>%
  left_join(tax.data) %>%
  mutate(hgcA = (binID %in% hgcA.list)) %>%
  left_join(elevation.info)



####-------------Set up color vector-------------####
color.vector <- cb.translator[c("vermillion", "yellow", "bluishgreen")]
names(color.vector) <- c("Kiritimatiellae-high", "Kiritimatiellae-low", "Lentisphaeria-low")
shape.vector <- c(15, 4)
names(shape.vector) <- c(TRUE, FALSE)


####-------------Plot out profiles of all bins-------------####
all.data %>%
  filter(((RM %in% c(286, 300)) & (year(date) %in% c(2017, 2018))) |
           ((RM %in% c(300, 310)) & (year(date) == 2019)))%>%
  ggplot(aes(y = coverage,
             x = elevation_m,
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
  coord_flip(xlim = c(545, 632)) +
  theme_bw()



####-------------Plot profiles of hgcA+ bins-------------####
all.data <- all.data %>%
  mutate(plotting.label = binID)
all.data[!all.data$hgcA, "plotting.label"] <- paste("hgcA- ",
                                                    all.data[!all.data$hgcA, "class"],
                                                    sep = "")


color.vector <- c(cb.translator[c("vermillion", "orange")],
                  "gray80", "gray20")
names(color.vector) <- c('anvio_hgcA_0261', 'anvio_hgcA_0040',
                         'hgcA- Kiritimatiellae')
KIR.2017.plots <- all.data %>%
  filter(RM %in% c(286, 300),
         year(date) == 2017,
         class == 'Kiritimatiellae') %>%
  arrange(desc(binID)) %>%
  ggplot(aes(y = coverage,
             x = elevation_m,
             group = binID,
             col = plotting.label)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = color.vector,
                     name = "") +
  scale_linetype_manual() +
  # scale_shape_manual(values = shape.vector,
  #                    name = "hgcA present") +
  scale_y_continuous(trans = "log10",
                     limits = c(0.0001, 3.5)) +
  coord_flip(xlim = c(545, 632)) +
  facet_wrap(~RM,
             nrow = 1, ncol = 2) +
  xlab("Depth (m)") +
  ylab("Bin abundance (%)") +
  theme_classic()
pdf("results/bins/binAnalysis/PVC_details/abundance_KIR_2017.pdf",
    height = 2.5,
    width = 4)
KIR.2017.plots
dev.off()  


####-------------2018 Kiritimatiellaeota-------------####
color.vector <- c(cb.translator[c("vermillion", "orange")], "gray80")
names(color.vector) <- c('HC18HY300_bin_0028', 'anvio_hgcA_0110', 'hgcA- Kiritimatiellae')
KIR.2018.plots <- all.data %>%
  filter(RM %in% c(286, 300),
         year(date) == 2018,
         class == "Kiritimatiellae") %>%
  arrange(desc(binID)) %>%
  ggplot(aes(y = coverage,
             x = elevation_m,
             group = binID,
             col = plotting.label)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = color.vector,
                     name = "") +
  scale_linetype_manual() +
  # scale_shape_manual(values = shape.vector,
  #                    name = "hgcA present") +
  scale_y_continuous(trans = "log10",
                     limits = c(0.0001, 3.5)) +
  coord_flip(xlim = c(545, 632)) +
  facet_wrap(~RM,
             nrow = 1, ncol = 2) +
  xlab("Elevation (m)") +
  ylab("Bin abundance (%)") +
  theme_classic()
pdf("results/bins/binAnalysis/PVC_details/abundance_KIR_2018.pdf",
    height = 2.5,
    width = 4)
KIR.2018.plots
dev.off()  



####-------------2017 Lentisphaerae-------------####
color.vector <- c(cb.translator["yellow"], "gray50")
names(color.vector) <- c('anvio_hgcA_0220', 'hgcA- Lentisphaeria')
LEN.2017.plots <- all.data %>%
  filter(RM %in% c(286, 300),
         year(date) == 2017,
         class == "Lentisphaeria") %>%
  arrange(desc(binID)) %>%
  ggplot(aes(y = coverage,
             x = elevation_m,
             group = binID,
             col = plotting.label)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = color.vector,
                     name = "") +
  scale_linetype_manual() +
  # scale_shape_manual(values = shape.vector,
  #                    name = "hgcA present") +
  scale_y_continuous(trans = "log10",
                     limits = c(0.0001, 3.5)) +
  coord_flip(xlim = c(545, 632)) +
  facet_wrap(~RM,
             nrow = 1, ncol = 2) +
  xlab("Elevation (m)") +
  ylab("Bin abundance (%)") +
  theme_classic()
pdf("results/bins/binAnalysis/PVC_details/abundance_LEN.2017.pdf",
    height = 2.5,
    width = 4)
LEN.2017.plots
dev.off()  
