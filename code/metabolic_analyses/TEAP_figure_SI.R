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


#### Read in depth data ####
gene.data <- readRDS("dataEdited/metabolic_analyses/depth/metabolicProtein_depth_clean.rds")
geochem.data <- rbind(read.csv("dataEdited/waterChemistry/geochem_WC.csv") %>%
                        mutate(date = as.Date(date)),
                      readRDS("dataEdited/seabird/seabird_data.rds"))


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
geochem.2017.286 <- plot.redox.profile(geochem.data.of.interest = geochem.data,
                                       RM.of.interest = 286,
                                       date.of.interest = "2017-09-25",
                                       depth.range.of.interest = c(75, 0),
                                       concentration.of.interest = c(0, 2),
                                       legend.location.of.interest = c(0.8, 0.5))
geochem.2017.300 <- plot.redox.profile(geochem.data.of.interest = geochem.data,
                                       RM.of.interest = 300,
                                       date.of.interest = "2017-09-26",
                                       depth.range.of.interest = c(75, 0),
                                       concentration.of.interest = c(0, 2),
                                       legend.location.of.interest = "none")
geochem.2018.286 <- plot.redox.profile(geochem.data.of.interest = geochem.data,
                                       RM.of.interest = 286,
                                       date.of.interest = "2018-09-24",
                                       depth.range.of.interest = c(75, 0),
                                       concentration.of.interest = c(0, 2),
                                       legend.location.of.interest = "none")
geochem.2018.300 <- plot.redox.profile(geochem.data.of.interest = geochem.data,
                                       RM.of.interest = 300,
                                       date.of.interest = "2018-09-25",
                                       depth.range.of.interest = c(75, 0),
                                       concentration.of.interest = c(0, 2),
                                       legend.location.of.interest = "none")
geochem.2019.300 <- plot.redox.profile(geochem.data.of.interest = geochem.data,
                                       RM.of.interest = 300,
                                       date.of.interest = "2019-07-25",
                                       depth.range.of.interest = c(60, 0),
                                       concentration.of.interest = c(0, 2),
                                       legend.location.of.interest = "none")
geochem.2019.310 <- plot.redox.profile(geochem.data.of.interest = geochem.data,
                                       RM.of.interest = 310,
                                       date.of.interest = "2019-07-23",
                                       depth.range.of.interest = c(60, 0),
                                       concentration.of.interest = c(0, 2),
                                       legend.location.of.interest = "none")
# ggarrange(geochem.2017.286, geochem.2017.300,
#           geochem.2018.286, geochem.2018.300,
#           geochem.2019.300, geochem.2019.310,
#           ncol = 2, nrow = 3)


#### Nitrate reduction ####
color.vector <- c(cb.translator["vermillion"],
                  cb.translator["orange"],
                  cb.translator["reddishpurple"])
names(color.vector) <- c("narG", "nirK", "nirS")
N.max.coverage <- 125
nar.286.2017 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                                genesOfInterest = names(color.vector),
                                                yearOfInterest = "2017",
                                                RMofInterest = "286",
                                                show.mean.coverage = FALSE,
                                                depth_limits = c(75, 0),
                                                coverage_limits = c(0, N.max.coverage),
                                                color.vector.to.use = color.vector,
                                                xlab.to.use = "Gene abundance",
                                                legend.position.to.use = c(0.8, 0.8))

nar.300.2017 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                                genesOfInterest = names(color.vector),
                                                yearOfInterest = "2017",
                                                RMofInterest = "300",
                                                show.mean.coverage = FALSE,
                                                depth_limits = c(75, 0),
                                                coverage_limits = c(0, N.max.coverage),
                                                color.vector.to.use = color.vector,
                                                xlab.to.use = "Gene abundance",
                                                legend.position.to.use = c(0.8, 0.8))

nar.286.2018 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                                genesOfInterest = names(color.vector),
                                                yearOfInterest = "2018",
                                                RMofInterest = "286",
                                                show.mean.coverage = FALSE,
                                                depth_limits = c(75, 0),
                                                coverage_limits = c(0, N.max.coverage),
                                                color.vector.to.use = color.vector,
                                                xlab.to.use = "Gene abundance",
                                                legend.position.to.use = "none")

nar.300.2018 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                                genesOfInterest = names(color.vector),
                                                yearOfInterest = "2018",
                                                RMofInterest = "300",
                                                show.mean.coverage = FALSE,
                                                depth_limits = c(75, 0),
                                                coverage_limits = c(0, N.max.coverage),
                                                color.vector.to.use = color.vector,
                                                xlab.to.use = "Gene abundance",
                                                legend.position.to.use = "none")
nar.300.2019 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                                genesOfInterest = names(color.vector),
                                                yearOfInterest = "2019",
                                                RMofInterest = "300",
                                                show.mean.coverage = FALSE,
                                                depth_limits = c(60, 0),
                                                coverage_limits = c(0, N.max.coverage),
                                                color.vector.to.use = color.vector,
                                                xlab.to.use = "Gene abundance",
                                                legend.position.to.use = "none")
nar.310.2019 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                                genesOfInterest = names(color.vector),
                                                yearOfInterest = "2019",
                                                RMofInterest = "310",
                                                show.mean.coverage = FALSE,
                                                depth_limits = c(60, 0),
                                                coverage_limits = c(0, N.max.coverage),
                                                color.vector.to.use = color.vector,
                                                xlab.to.use = "Gene abundance",
                                                legend.position.to.use = "none")




#### Sulfate cycling ####
color.vector <- c(cb.translator["blue"], cb.translator["skyblue"])
names(color.vector) <- c("dsrA", "mcrA")
max.coverage <- 4.5
S.286.2017 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                              genesOfInterest = names(color.vector),
                                              yearOfInterest = "2017",
                                              RMofInterest = "286",
                                              show.mean.coverage = FALSE,
                                              depth_limits = c(75, 0),
                                              coverage_limits = c(0, max.coverage),
                                              color.vector.to.use = color.vector,
                                              xlab.to.use = "Gene abundance",
                                              legend.position.to.use = c(0.8, 0.8))

S.300.2017 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                              genesOfInterest = names(color.vector),
                                              yearOfInterest = "2017",
                                              RMofInterest = "300",
                                              show.mean.coverage = FALSE,
                                              depth_limits = c(75, 0),
                                              coverage_limits = c(0, max.coverage),
                                              color.vector.to.use = color.vector,
                                              xlab.to.use = "Gene abundance",
                                              legend.position.to.use = c(0.8, 0.8))

S.286.2018 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                              genesOfInterest = names(color.vector),
                                              yearOfInterest = "2018",
                                              RMofInterest = "286",
                                              show.mean.coverage = FALSE,
                                              depth_limits = c(75, 0),
                                              coverage_limits = c(0, max.coverage),
                                              color.vector.to.use = color.vector,
                                              xlab.to.use = "Gene abundance",
                                              legend.position.to.use = "none")

S.300.2018 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                              genesOfInterest = names(color.vector),
                                              yearOfInterest = "2018",
                                              RMofInterest = "300",
                                              show.mean.coverage = FALSE,
                                              depth_limits = c(75, 0),
                                              coverage_limits = c(0, max.coverage),
                                              color.vector.to.use = color.vector,
                                              xlab.to.use = "Gene abundance",
                                              legend.position.to.use = "none")

S.300.2019 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                              genesOfInterest = names(color.vector),
                                              yearOfInterest = "2019",
                                              RMofInterest = "300",
                                              show.mean.coverage = FALSE,
                                              depth_limits = c(60, 0),
                                              coverage_limits = c(0, max.coverage),
                                              color.vector.to.use = color.vector,
                                              xlab.to.use = "Gene abundance",
                                              legend.position.to.use = "none")

S.310.2019 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                              genesOfInterest = names(color.vector),
                                              yearOfInterest = "2019",
                                              RMofInterest = "310",
                                              show.mean.coverage = FALSE,
                                              depth_limits = c(60, 0),
                                              coverage_limits = c(0, max.coverage),
                                              color.vector.to.use = color.vector,
                                              xlab.to.use = "Gene abundance",
                                              legend.position.to.use = "none")







pdf("results/manuscript_figures/TEAP_SI_figure.pdf",
    width = 7.0,
    height = 8)

ggarrange(geochem.2017.286, nar.286.2017, S.286.2017,
          geochem.2017.300, nar.300.2017, S.300.2017,
          geochem.2018.286, nar.286.2018, S.286.2018,
          geochem.2018.300, nar.300.2018, S.300.2018,
          geochem.2019.300, nar.300.2019, S.300.2019,
          geochem.2019.310, nar.310.2019, S.310.2019,
          ncol = 6, nrow = 3
          )
dev.off()
