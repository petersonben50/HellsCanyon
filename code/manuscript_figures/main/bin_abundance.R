#### code/manuscript_figures/main/bin_abundance.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
cb.translator <- c(cb.translator, "gray80")
names(cb.translator)[length(cb.translator)] <- "gray80"
cb.translator <- c(cb.translator, "gray50")
names(cb.translator)[length(cb.translator)] <- "gray50"
cb.translator <- c(cb.translator, "gray20")
names(cb.translator)[length(cb.translator)] <- "gray20"


#### Read in metadata ####
MG.metadata <- read_xlsx("metadata/metagenome_metadata.xlsx")


#### Read in metabolism data ####
metabolic.data <- read_xlsx("dataEdited/binning/metabolism/metabolic_summary.xlsx",
                            sheet = "overall_summary") %>%
  select(binID, assigned_group)
HMS.plotting <- read_xlsx("dataEdited/binning/metabolism/metabolic_summary.xlsx",
                        sheet = "HMS_plotting")
HMS.colors.vector <- cb.translator[HMS.plotting$colorToUse]
names(HMS.colors.vector) <- HMS.plotting$binID
HMS.points.vector <- HMS.plotting$pointToUse
names(HMS.points.vector) <- HMS.plotting$binID


#### Read in geochem data ####
geochem.data <- read.csv("dataEdited/waterChemistry/geochem_WC.csv") %>%
  group_by(RM, depth, date, constituent) %>%
  summarize(concentration = mean(concentration)) %>%
  spread(key = constituent,
         value = concentration)



#### Read in depth data ####
depth.data <- readRDS("dataEdited/binning/depth/bin_depth_clean.rds") %>%
  left_join(MG.metadata %>% select(metagenomeID, date, RM, depth)) %>%
  left_join(metabolic.data) %>%
  filter(!is.na(assigned_group))
rm(metabolic.data, MG.metadata)


#### Generate plot for fermenters ####
pdf("results/manuscript_figures/binning_fig/fermenters_2017.pdf",
    width = 2.25,
    height = 2.25)
depth.data %>%
  filter(assigned_group == "fermenter") %>%
  filter(year(date) == 2017) %>%
  ggplot(aes(y = coverage,
             x = depth,
             group = binID)) +
  geom_point(aes(color = binID,
                 pch = binID)) +
  geom_line(aes(color = binID)) +
  scale_color_manual(values = HMS.colors.vector) +
  scale_shape_manual(values = HMS.points.vector) +
  facet_wrap(~RM, nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank()) +
  labs(y = "Gene abundance",
       x = "Depth (m)",
       title = element_blank()) +
  coord_flip(xlim = c(80, 0))
dev.off()

pdf("results/manuscript_figures/binning_fig/fermenters_2018.pdf",
    width = 2.25,
    height = 2.25)
depth.data %>%
  filter(assigned_group == "fermenter") %>%
  filter(year(date) == 2018) %>%
  ggplot(aes(y = coverage,
             x = depth,
             group = binID)) +
  geom_point(aes(color = binID,
                 pch = binID)) +
  geom_line(aes(color = binID)) +
  scale_color_manual(values = HMS.colors.vector) +
  scale_shape_manual(values = HMS.points.vector) +
  facet_wrap(~RM, nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank()) +
  labs(y = "Gene abundance",
       x = "Depth (m)",
       title = element_blank()) +
  coord_flip(xlim = c(80, 0))
dev.off()




#### Generate plots for HRR ####
pdf("results/manuscript_figures/binning_fig/HRR_2018.pdf",
    width = 2.25,
    height = 2.25)
depth.data %>%
  filter(assigned_group == "HRR") %>%
  filter(year(date) == 2018) %>%
  ggplot(aes(y = coverage,
             x = depth,
             group = binID)) +
  geom_point(aes(color = binID,
                 pch = binID)) +
  geom_line(aes(color = binID)) +
  scale_color_manual(values = HMS.colors.vector) +
  scale_shape_manual(values = HMS.points.vector) +
  facet_wrap(~RM, nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank()) +
  labs(y = "Gene abundance",
       x = "Depth (m)",
       title = element_blank()) +
  coord_flip(xlim = c(80, 0))
dev.off()

pdf("results/manuscript_figures/binning_fig/HRR_2019.pdf",
    width = 2.25,
    height = 2.25)
depth.data %>%
  filter(assigned_group == "HRR") %>%
  filter(year(date) == 2019) %>%
  filter(RM %in% c(300, 310)) %>%
  ggplot(aes(y = coverage,
             x = depth,
             group = binID)) +
  geom_point(aes(color = binID,
                 pch = binID)) +
  geom_line(aes(color = binID)) +
  scale_color_manual(values = HMS.colors.vector) +
  scale_shape_manual(values = HMS.points.vector) +
  facet_wrap(~RM, nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank()) +
  labs(y = "Gene abundance",
       x = "Depth (m)",
       title = element_blank()) +
  coord_flip(xlim = c(80, 0))
dev.off()




#### Generate plots for SRB ####
pdf("results/manuscript_figures/binning_fig/SRB_2017.pdf",
    width = 2.25,
    height = 2.25)
depth.data %>%
  filter(assigned_group == "SRB") %>%
  filter(year(date) == 2017) %>%
  ggplot(aes(y = coverage,
             x = depth,
             group = binID)) +
  geom_point(aes(color = binID,
                 pch = binID)) +
  geom_line(aes(color = binID)) +
  scale_color_manual(values = HMS.colors.vector) +
  scale_shape_manual(values = HMS.points.vector) +
  facet_wrap(~RM, nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank()) +
  labs(y = "Gene abundance",
       x = "Depth (m)",
       title = element_blank()) +
  coord_flip(xlim = c(80, 0))
dev.off()

pdf("results/manuscript_figures/binning_fig/SRB_2019.pdf",
    width = 2.25,
    height = 2.25)
depth.data %>%
  filter(assigned_group == "SRB") %>%
  filter(year(date) == 2019) %>%
  filter(RM %in% c(300, 310)) %>%
  ggplot(aes(y = coverage,
             x = depth,
             group = binID)) +
  geom_point(aes(color = binID,
                 pch = binID)) +
  geom_line(aes(color = binID)) +
  scale_color_manual(values = HMS.colors.vector) +
  scale_shape_manual(values = HMS.points.vector) +
  facet_wrap(~RM, nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank()) +
  labs(y = "Gene abundance",
       x = "Depth (m)",
       title = element_blank()) +
  coord_flip(xlim = c(80, 0))
dev.off()