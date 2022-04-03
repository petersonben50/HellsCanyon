#### code/metabolic_analyses/TEAP_figure.R ####
# Benjamin D. Peterson

# This script generates the figure that plots
# DO, sulfide, nitrate, and Mn for comparison
# to denitrification, sulfate reduction, and 
# methanogenesis genes. 

#### Clean up crew on line 5 ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(ggpubr)
library(lubridate)
library(patchwork)
library(tidyverse)
source("code/plotting_functions.R")
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Read in data ####
gene.data <- readRDS("dataEdited/metabolic_analyses/depth/metabolicProtein_depth_clean.rds")
geochem.data <- rbind(read.csv("dataEdited/waterChemistry/geochem_WC.csv") %>%
                        mutate(date = as.Date(date)),
                      readRDS("dataEdited/seabird/seabird_data.rds"))
hgcA.data <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv") %>%
  filter(rep) %>%
  group_by(date, RM, depth) %>%
  summarise(coverage = sum(coverage))
# Set 0.01 as the lower limit
hgcA.data[hgcA.data$coverage < 0.01, "coverage"] <- 0.011

#### Vectors for geochem ####
color.vector <- c(cb.translator["black"],
                  cb.translator["vermillion"],
                  cb.translator["orange"],
                  cb.translator["blue"])
names(color.vector) <- c("DO", "f_no3_mg_n_per_l", "f_mn_mg_per_l", "f_inorganic_sulfide_mg_per_l")
line.vector <- c(4, 1, 3, 2)
names(line.vector) <- names(color.vector)
labels.vector <- c("DO (mg/L)", "NO3 (mgN/L)", "Diss. Mn", "Sulfide")
names(labels.vector) <- names(color.vector)
points.vector <- c(NA, 21, 24, 23)
names(points.vector) <- names(color.vector)
scaling.vector <- 0.25
names(scaling.vector) <- "DO"



#### Generate geochem profiles ####
geochem.2017.286 <- plot.DO.nitrate.Mn.sulfide.profile(geochem.data.of.interest = geochem.data,
                                                       RM.of.interest = 286,
                                                       date.of.interest = "2017-09-25",
                                                       depth.range.of.interest = c(75, 0),
                                                       concentration.of.interest = c(0, 2),
                                                       legend.location.of.interest = c(0.8, 0.5),
                                                       scaling.vector.to.use = scaling.vector)
geochem.2017.300 <- plot.DO.nitrate.Mn.sulfide.profile(geochem.data.of.interest = geochem.data,
                                                       RM.of.interest = 300,
                                                       date.of.interest = "2017-09-26",
                                                       depth.range.of.interest = c(75, 0),
                                                       concentration.of.interest = c(0, 2),
                                                       legend.location.of.interest = "none",
                                                       scaling.vector.to.use = scaling.vector)
geochem.2018.286 <- plot.DO.nitrate.Mn.sulfide.profile(geochem.data.of.interest = geochem.data,
                                                       RM.of.interest = 286,
                                                       date.of.interest = "2018-09-24",
                                                       depth.range.of.interest = c(75, 0),
                                                       concentration.of.interest = c(0, 2),
                                                       legend.location.of.interest = "none",
                                                       scaling.vector.to.use = scaling.vector)
geochem.2018.300 <- plot.DO.nitrate.Mn.sulfide.profile(geochem.data.of.interest = geochem.data,
                                                       RM.of.interest = 300,
                                                       date.of.interest = "2018-09-25",
                                                       depth.range.of.interest = c(75, 0),
                                                       concentration.of.interest = c(0, 2),
                                                       legend.location.of.interest = "none",
                                                       scaling.vector.to.use = scaling.vector)



#### TEAP genes ####
color.vector <- c(cb.translator["vermillion"],
                  cb.translator["blue"],
                  cb.translator["reddishpurple"]
                  )
names(color.vector) <- c("narG", "dsrA", "mcrA")
points.vector <- c(21, 23, 22)
names(points.vector) <- names(color.vector)
lines.vector <- c(1, 2, 6)
names(lines.vector) <- names(color.vector)

max.coverage <- 50
TEAPs.286.2017 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                                  genesOfInterest = names(color.vector),
                                                  yearOfInterest = "2017",
                                                  RMofInterest = "286",
                                                  show.mean.coverage = FALSE,
                                                  depth_limits = c(75, 0),
                                                  coverage_limits = c(0.01, max.coverage),
                                                  color.vector.to.use = color.vector,
                                                  xlab.to.use = "Gene abundance",
                                                  point.vector.to.use = points.vector,
                                                  line.vector.to.use = lines.vector,
                                                  DL = 0.01,
                                                  legend.position.to.use = c(0.8, 0.8))

TEAPs.300.2017 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                                  genesOfInterest = names(color.vector),
                                                  yearOfInterest = "2017",
                                                  RMofInterest = "300",
                                                  show.mean.coverage = FALSE,
                                                  depth_limits = c(75, 0),
                                                  coverage_limits = c(0.01, max.coverage),
                                                  color.vector.to.use = color.vector,
                                                  xlab.to.use = "Gene abundance",
                                                  point.vector.to.use = points.vector,
                                                  line.vector.to.use = lines.vector,
                                                  DL = 0.01,
                                                  legend.position.to.use = c(0.8, 0.8))

TEAPs.286.2018 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                                  genesOfInterest = names(color.vector),
                                                  yearOfInterest = "2018",
                                                  RMofInterest = "286",
                                                  show.mean.coverage = FALSE,
                                                  depth_limits = c(75, 0),
                                                  coverage_limits = c(0.01, max.coverage),
                                                  color.vector.to.use = color.vector,
                                                  xlab.to.use = "Gene abundance",
                                                  point.vector.to.use = points.vector,
                                                  line.vector.to.use = lines.vector,
                                                  DL = 0.01,
                                                  legend.position.to.use = "none")

TEAPs.300.2018 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                                  genesOfInterest = names(color.vector),
                                                  yearOfInterest = "2018",
                                                  RMofInterest = "300",
                                                  show.mean.coverage = FALSE,
                                                  depth_limits = c(75, 0),
                                                  coverage_limits = c(0.01, max.coverage),
                                                  color.vector.to.use = color.vector,
                                                  xlab.to.use = "Gene abundance",
                                                  point.vector.to.use = points.vector,
                                                  line.vector.to.use = lines.vector,
                                                  DL = 0.01,
                                                  legend.position.to.use = "none")




#### hgcA and MeHg genes ####
date.of.interest <- "2017-09-25"
RM.of.interest <- 286

color.vector <- c(cb.translator["skyblue"], cb.translator["black"])
names(color.vector) <- c("hgcA", "MeHg_diss_ngL")
points.vector <- c(21, 23)
names(points.vector) <- names(color.vector)
line.vector <- c(1, 2)
names(points.vector) <- names(line.vector)

Hg.data.profile <- rbind(hgcA.data %>%
                           rename(amount = coverage) %>%
                           mutate(constituent = "hgcA",
                                  date = as.Date(date),
                                  RM = as.character(RM)) %>%
                           select(date, RM, depth, constituent, amount),
                         geochem.data %>%
                           filter(constituent == "MeHg_diss_ngL") %>%
                           rename(amount = concentration))

Hg.286.2017 <- plot.hgcA.and.MeHg(data.to.use = Hg.data.profile,
                                  date.of.interest = "2017-09-25",
                                  RM.of.interest = 286,
                                  depth.range.of.interest = c(75, 0),
                                  legend.location.of.interest = c(0.8, 0.8),
                                  abundance.limits = c(0.005, 10))
Hg.300.2017 <- plot.hgcA.and.MeHg(data.to.use = Hg.data.profile,
                                  date.of.interest = "2017-09-26",
                                  RM.of.interest = 300,
                                  depth.range.of.interest = c(75, 0),
                                  legend.location.of.interest = c(0.8, 0.8),
                                  abundance.limits = c(0.005, 10))
Hg.286.2018 <- plot.hgcA.and.MeHg(data.to.use = Hg.data.profile,
                                  date.of.interest = "2018-09-24",
                                  RM.of.interest = 286,
                                  depth.range.of.interest = c(75, 0),
                                  legend.location.of.interest = "none",
                                  abundance.limits = c(0.005, 10))
Hg.300.2018 <- plot.hgcA.and.MeHg(data.to.use = Hg.data.profile,
                                  date.of.interest = "2018-09-25",
                                  RM.of.interest = 300,
                                  depth.range.of.interest = c(75, 0),
                                  legend.location.of.interest = "none",
                                  abundance.limits = c(0.005, 10))



pdf("results/manuscript_figures/BGC_main_figure.pdf",
    width = 11,
    height = 8)

ggarrange(geochem.2017.286, TEAPs.286.2017, Hg.286.2017,
          geochem.2017.300, TEAPs.300.2017, Hg.300.2017,
          geochem.2018.286, TEAPs.286.2018, Hg.286.2018,
          geochem.2018.300, TEAPs.300.2018, Hg.300.2018,
          ncol = 6, nrow = 2
          )
dev.off()
