#### code/binning/bin_depth/SRB_abundance.R ####
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
MG.metadata <- read_xlsx("metadata/metagenome_metadata.xlsx") %>%
  left_join(read.csv("dataEdited/geochem/geochem_WC.csv") %>%
              select(date, RM, depth, elevation_m) %>%
              mutate(RM = as.character(RM)) %>%
              unique())

#### Read in metabolism data ####
metabolic.data <- read_xlsx("dataEdited/bins/binAnalysis/metabolism/aggregated_metabolism.xlsx",
                            sheet = "Metabolism_summary_hgcA_bins") %>%
  rename(metabolic_assignment = `Metabolic classification`) %>%
  rename(HMS = `mOTU ID`) %>%
  select(HMS, binID, metabolic_assignment)



#### Read in depth data ####
depth.data <- readRDS("dataEdited/bins/binAnalysis/depth/bin_depth_clean.rds") %>%
  left_join(MG.metadata %>% select(metagenomeID, date, RM, depth, elevation_m)) %>%
  left_join(metabolic.data)%>%
  filter(((RM %in% c(286, 300)) & (year(date) %in% c(2017, 2018))) |
           ((RM %in% c(300, 310)) & (year(date) == 2019))) %>%
  group_by(binID, RM, date, depth, elevation_m, metabolic_assignment) %>%
  summarise(coverage = mean(coverage))
rm(metabolic.data, MG.metadata)


#### Set up color vector ####
color.vector <- cb.translator[c("skyblue", "blue", "bluishgreen")]
names(color.vector) <- c("anvio_hgcA_0240", "anvio_hgcA_0070", "anvio_hgcA_0250")


#### Plotting function ####
plot.bins.by.year <- function(depth.data.to.use = depth.data,
                              year.of.interest,
                              color.vector.to.use = color.vector) {
  depth.data.to.use %>%
    filter(year(date) == year.of.interest,
           binID %in% names(color.vector.to.use)) %>%
    ggplot(aes(y = coverage,
             x = elevation_m,
             group = binID)) +
    geom_point(aes(color = binID)) +
    geom_line(aes(color = binID)) +
    scale_color_manual(values = color.vector.to.use) +
    facet_wrap(~RM, nrow = 1) +
    scale_y_continuous(limits = c(0, 1.5)) +
    coord_flip(xlim = c(545, 632)) +
    xlab("Elevation (m)") +
    ylab("Bin abundance (%)") +
    theme_classic()
}


pdf("results/bins/binAnalysis/depth/SRB_abundance.pdf",
    height = 2.5,
    width = 4)
plot.bins.by.year(year.of.interest = 2017)
dev.off()
plot.bins.by.year(year.of.interest = 2019)
