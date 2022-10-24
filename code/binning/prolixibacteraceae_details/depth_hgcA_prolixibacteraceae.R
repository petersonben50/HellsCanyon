#### code/binning/prolixibacteraceae_details/depth_hgcA_prolixibacteraceae.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(ggpubr)
library(lubridate)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")
cb.translator <- c(cb.translator, "grey50")
names(cb.translator)[length(cb.translator)] <- "grey"


#### Read in metadata ####
MG.metadata <- read_xlsx("metadata/metagenome_metadata.xlsx")



#### Read in depth data ####
depth.data <- readRDS("dataEdited/bins/binAnalysis/depth/bin_depth_clean.rds") %>%
  left_join(MG.metadata %>% select(metagenomeID, date, RM, depth)) %>%
  filter(binID == "anvio_hgcA_0130")



#### Generate depth plots ####
bin.depth.286 <- depth.data %>%
  filter(RM == "286") %>%
  ggplot(aes(y = coverage,
             x = depth)) +
  geom_point(color = cb.translator["orange"]) +
  geom_line(color = cb.translator["orange"]) +
  labs(y = "Genome\nabundance (%)",
       x = "Depth (m)",
       title = "Sept 2018 - RM286") +
  coord_flip(xlim = c(80, 0),
             ylim = c(0, 0.3)) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

bin.depth.300 <- depth.data %>%
  filter(RM == "300") %>%
  ggplot(aes(y = coverage,
             x = depth)) +
  geom_point(color = cb.translator["orange"]) +
  geom_line(color = cb.translator["orange"]) +
  labs(y = "Genome\nabundance (%)",
       x = "Depth (m)",
       title = "Sept 2018 - RM300") +
  coord_flip(xlim = c(80, 0),
             ylim = c(0, 0.3)) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))


#### Save out plots ####
pdf("results/bins/binAnalysis/prolixibacteraceae_details/depth_anvio_hgcA_0130.pdf",
    height = 3,
    width = 3)
ggarrange(bin.depth.286,
          bin.depth.300)
dev.off()
