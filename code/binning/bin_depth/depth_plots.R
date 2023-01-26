#### code/binning/depth_plots.R ####
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


#### Read in metabolism data ####
metabolic.data <- read_xlsx("dataEdited/bins/binAnalysis/metabolism/aggregated_metabolism.xlsx",
                            sheet = "Metabolism_summary_hgcA_bins") %>%
  rename(metabolic_assignment = `Metabolic classification`) %>%
  rename(HMS = `mOTU ID`) %>%
  select(HMS, binID, metabolic_assignment)



#### Prepare colors ####
metabolism.colors <- read_xlsx("dataEdited/bins/binAnalysis/metabolism/aggregated_metabolism.xlsx",
                               sheet = "colors_to_use")
metabolism.colors.vector <- cb.translator[metabolism.colors$color_to_use]
names(metabolism.colors.vector) <- metabolism.colors$metabolic_assignment
rm(metabolism.colors, cb.translator)


#### Read in depth data ####
depth.data <- readRDS("dataEdited/bins/binAnalysis/depth/bin_depth_clean.rds") %>%
  left_join(MG.metadata %>% select(metagenomeID, date, RM, depth)) %>%
  left_join(metabolic.data)%>%
  filter(((RM %in% c(286, 300)) & (year(date) %in% c(2017, 2018))) |
           ((RM %in% c(300, 310)) & (year(date) == 2019))) %>%
  group_by(HMS, RM, date, depth, metabolic_assignment) %>%
  summarise(coverage = mean(coverage))
rm(metabolic.data, MG.metadata)


#### Plotting function ####
plot.bins.by.year <- function(year.of.interest,
                              metabolism_of_interest = NULL) {
  if (!is.null(metabolism_of_interest)) {
    depth.data.to.use <- depth.data %>%
      filter(metabolic_assignment %in% metabolism_of_interest)
    metabolism.colors.vector <- metabolism.colors.vector[metabolism_of_interest]
  } else {
    depth.data.to.use <- depth.data
  }
  depth.data.to.use %>%
    filter(year(date) == year.of.interest) %>%
    ggplot(aes(y = coverage,
             x = depth,
             group = HMS)) +
    geom_point(aes(color = metabolic_assignment)) +
    geom_line(aes(color = metabolic_assignment)) +
    scale_color_manual(values = metabolism.colors.vector) +
    facet_wrap(~RM, nrow = 1) +
    coord_flip(xlim = c(80, 0)) +
    theme_bw()
}


plot.bins.by.year(2017)
bins.2017 <- plot.bins.by.year(2017, metabolism_of_interest = c("fermentative", "fermentative_O2_tolerant"))
bins.2018 <- plot.bins.by.year(2018, metabolism_of_interest = c("fermentative", "fermentative_O2_tolerant"))
ggarrange(bins.2017, bins.2018,
          ncol = 1)


plot.bins.by.year(2017, metabolism_of_interest = c("sulfate-reducer"))
