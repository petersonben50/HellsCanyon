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
geochem.data <- read.csv("dataEdited/waterChemistry/geochem_WC.csv") %>%
                        mutate(date = as.Date(date))
seabird.data <- readRDS("dataEdited/seabird/seabird_data.rds")



#### Generate sonde profiles ####

# Vectors for sonde profiles
color.vector <- c(cb.translator["skyblue"],
                  cb.translator["black"])
names(color.vector) <- c("temp", "DO")
line.vector <- c(1, 3)
names(line.vector) <- c("temp", "DO")
labels.vector <- c("Temperature", "Dissolved O2 (mg/L) (2.5X)")
names(labels.vector) <- c("temp", "DO")
scaling.vector <- 2.5
names(scaling.vector) <- "DO"

sonde.2019.300 <- ggplot() + theme_void()
sonde.2019.310 <- sonde.profile(seabird.data.of.interest = seabird.data,
                                RM.of.interest = 310,
                                year.of.interest =  2019,
                                depth.range.of.interest = c(60, 0),
                                concentration.of.interest = c(0, 25),
                                legend.location.of.interest = "none",
                                scaling.vector.to.use = scaling.vector,
                                color.vector.to.use = color.vector)



#### Vectors for geochem ####
color.vector <- c(cb.translator["reddishpurple"],
                  cb.translator["orange"],
                  cb.translator["blue"])
names(color.vector) <- c("f_no3_mg_n_per_l", "f_mn_mg_per_l", "f_inorganic_sulfide_mg_per_l")
line.vector <- c(1, 3, 2)
names(line.vector) <- names(color.vector)
labels.vector <- c("NO3 (mgN/L)", "Diss. Mn", "Sulfide")
names(labels.vector) <- names(color.vector)
points.vector <- c(21, 24, 23)
names(points.vector) <- names(color.vector)


#### Generate geochem profiles ####
geochem.2019.300 <- plot.redox.profile(geochem.data.of.interest = geochem.data,
                                       RM.of.interest = 300,
                                       date.of.interest = "2019-07-25",
                                       depth.range.of.interest = c(60, 0),
                                       concentration.of.interest = c(0, 2),
                                       legend.location.of.interest = c(0.8, 0.5))
geochem.2019.310 <- plot.redox.profile(geochem.data.of.interest = geochem.data,
                                       RM.of.interest = 310,
                                       date.of.interest = "2019-07-23",
                                       depth.range.of.interest = c(60, 0),
                                       concentration.of.interest = c(0, 2),
                                       legend.location.of.interest = "none")



#### Hg speciation profiles ####
date.of.interest <- "2017-09-25"
RM.of.interest <- 286
color.vector <- c(cb.translator["bluishgreen"], cb.translator["vermillion"], cb.translator["black"])
names(color.vector) <- c("iHg_diss_ngL", "MeHg_diss_ngL", "MeHg_diss_percent")
points.vector <- c(21, 24, 23)
names(points.vector) <- names(color.vector)
labels.vector <- c("HgT (ng/L)", "MeHg (ng/L)", "MeHg (%) (0.05X")
names(labels.vector) <- names(color.vector)
line.vector <- c(1, 3, 2)
names(points.vector) <- names(line.vector)
scaling.vector <- 0.05
names(scaling.vector) <- "MeHg_diss_percent"

Hg.concentration.range.to.use <- c(0, 4)

Hg.2019.300 <- plot.redox.profile(geochem.data.of.interest = geochem.data,
                                  RM.of.interest = 300,
                                  date.of.interest = "2019-07-25",
                                  depth.range.of.interest = c(60, 0),
                                  concentration.of.interest = Hg.concentration.range.to.use,
                                  legend.location.of.interest = c(0.8, 0.5),
                                  scaling.vector.to.use = scaling.vector)
Hg.2019.310 <- plot.redox.profile(geochem.data.of.interest = geochem.data,
                                  RM.of.interest = 310,
                                  date.of.interest = "2019-07-23",
                                  depth.range.of.interest = c(60, 0),
                                  concentration.of.interest = Hg.concentration.range.to.use,
                                  legend.location.of.interest = "none",
                                  scaling.vector.to.use = scaling.vector)



#### Produce plot ####
pdf("results/manuscript_figures/BGC_2019_SI_figure.pdf",
    width = 11,
    height = 4)

ggarrange(sonde.2019.300, geochem.2019.300, Hg.2019.300,
          sonde.2019.310, geochem.2019.310, Hg.2019.310,
          ncol = 6, nrow = 1
          )
dev.off()
