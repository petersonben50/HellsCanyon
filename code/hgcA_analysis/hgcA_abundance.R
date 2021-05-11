#### code/hgcA_analysis/hgcA_abundance.R ####
# Written by Benjamin D. Peterson


#### Get ready. Get set ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(lubridate)
library(patchwork)
library(readxl)
library(tidyverse)
source("code/gene_plotting_functions.R")
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
cb.translator <- c(cb.translator, "gray50")
names(cb.translator)[length(cb.translator)] <- "gray50"


#### Read in hgcA classification ####
tax.data <- read_xlsx("dataEdited/hgcA_analysis/hgcA_information_edited.xlsx") %>%
  select(seqID, manual_classification, estimatedMetabolism)


#### Make color vector ####
hgcA.manual.taxonomy <- read_xlsx("dataEdited/hgcA_analysis/phylogeny/manual_taxonomy.xlsx",
                                  sheet = "colors_to_use")
color.vector <- cb.translator[hgcA.manual.taxonomy$colorsToUse]
names(color.vector) <- hgcA.manual.taxonomy$seqID


#### Read in data ####
hgcA.data <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv") %>%
  filter(rep) %>%
  left_join(tax.data)


#### 2017 profiles ####
profile.286.2017 <- plot.profile.of.gene.with.taxonomy(data.to.use = hgcA.data,
                                                       RM.of.interest = 286,
                                                       year.of.interest = 2017,
                                                       legend.position.to.use = "none",
                                                       coverage_limits = c(0, 10),
                                                       taxonomy.column.name = "manual_classification",
                                                       ylab.to.use = "hgcA gene coverage\n(per 100X SCG)",
                                                       titleToUse = "RM286 in 2017")
profile.300.2017 <- plot.profile.of.gene.with.taxonomy(data.to.use = hgcA.data,
                                                       RM.of.interest = 300,
                                                       year.of.interest = 2017,
                                                       legend.position.to.use = c(0.7, 0.8),
                                                       coverage_limits = c(0, 2.5),
                                                       taxonomy.column.name = "manual_classification",
                                                       ylab.to.use = "hgcA gene coverage\n(per 100X SCG)",
                                                       titleToUse = "RM300 in 2017")
pdf("results/hgcA_analysis/profiles/2017_profiles.pdf",
    height = 4.5,
    width = 6)
profile.286.2017 + profile.300.2017
dev.off()

#### 2017 close-up profile of metalimnion at 286 ####
profile.286.2017.closeup <- plot.profile.of.gene.with.taxonomy(data.to.use = hgcA.data,
                                                       RM.of.interest = 286,
                                                       year.of.interest = 2017,
                                                       legend.position.to.use = c(0.8, 0.8),
                                                       coverage_limits = c(0, 0.25),
                                                       depth_limits = c(40, 18),
                                                       taxonomy.column.name = "manual_classification",
                                                       ylab.to.use = "hgcA gene coverage")
pdf("results/hgcA_analysis/profiles/2017_profiles_286metalimnion.pdf",
    height = 4,
    width = 4)
profile.286.2017.closeup
dev.off()


#### 2017 close-up profile of metalimnion at 286, show inferred metabolism ####
color.vector.metabolism <- c(cb.translator["vermillion"], cb.translator["bluishgreen"], cb.translator["orange"], cb.translator["reddishpurple"], cb.translator["black"])
names(color.vector.metabolism) <- c("Fermenter", "SRB", "MRB", "methanogen", "unknown")

profile.286.2017.closeup <- plot.profile.of.gene.with.taxonomy(data.to.use = hgcA.data,
                                                               RM.of.interest = 286,
                                                               year.of.interest = 2017,
                                                               legend.position.to.use = c(0.8, 0.8),
                                                               color.vector.to.use = color.vector.metabolism,
                                                               coverage_limits = c(0, 0.25),
                                                               depth_limits = c(40, 18),
                                                               taxonomy.column.name = "estimatedMetabolism",
                                                               ylab.to.use = "hgcA gene coverage")
pdf("results/hgcA_analysis/profiles/2017_profiles_286metalimnion_metabolism.pdf",
    height = 4,
    width = 4)
profile.286.2017.closeup
dev.off()

#### 2017 close-up profile of hypolimnion at 286, show inferred metabolism ####
profile.286.2017.closeup <- plot.profile.of.gene.with.taxonomy(data.to.use = hgcA.data,
                                                               RM.of.interest = 286,
                                                               year.of.interest = 2017,
                                                               legend.position.to.use = c(0.8, 0.8),
                                                               color.vector.to.use = color.vector.metabolism,
                                                               coverage_limits = c(0, 10),
                                                               depth_limits = c(80, 40),
                                                               taxonomy.column.name = "estimatedMetabolism",
                                                               ylab.to.use = "hgcA gene coverage")
pdf("results/hgcA_analysis/profiles/2017_profiles_286hypolimnion_metabolism.pdf",
    height = 4,
    width = 4)
profile.286.2017.closeup
dev.off()

#### 2017 close-up profile of hypolimnion at 300, show inferred metabolism ####
profile.300.2017.closeup <- plot.profile.of.gene.with.taxonomy(data.to.use = hgcA.data,
                                                               RM.of.interest = 300,
                                                               year.of.interest = 2017,
                                                               legend.position.to.use = c(0.8, 0.8),
                                                               color.vector.to.use = color.vector.metabolism,
                                                               coverage_limits = c(0, 2),
                                                               depth_limits = c(55, 35),
                                                               taxonomy.column.name = "estimatedMetabolism",
                                                               ylab.to.use = "hgcA gene coverage")
pdf("results/hgcA_analysis/profiles/2017_profiles_300hypolimnion_metabolism.pdf",
    height = 4,
    width = 4)
profile.300.2017.closeup
dev.off()


#### 2018 profiles ####
profile.286.2018 <- plot.profile.of.gene.with.taxonomy(data.to.use = hgcA.data,
                                                       RM.of.interest = 286,
                                                       year.of.interest = 2018,
                                                       legend.position.to.use = "none",
                                                       coverage_limits = c(0, 1),
                                                       taxonomy.column.name = "manual_classification",
                                                       ylab.to.use = "hgcA gene coverage\n(per 100X SCG)",
                                                       titleToUse = "RM286 in 2018")
profile.300.2018 <- plot.profile.of.gene.with.taxonomy(data.to.use = hgcA.data,
                                                       RM.of.interest = 300,
                                                       year.of.interest = 2018,
                                                       legend.position.to.use = c(0.6, 0.8),
                                                       coverage_limits = c(0, 3.5),
                                                       taxonomy.column.name = "manual_classification",
                                                       ylab.to.use = "hgcA gene coverage\n(per 100X SCG)",
                                                       titleToUse = "RM300 in 2018")

pdf("results/hgcA_analysis/profiles/2018_profiles.pdf",
    height = 4.5,
    width = 6)
profile.286.2018 + profile.300.2018
dev.off()


#### 2018 close-up profile of metalimnion at 286, show inferred metabolism ####
profile.286.2018.closeup <- plot.profile.of.gene.with.taxonomy(data.to.use = hgcA.data,
                                                               RM.of.interest = 286,
                                                               year.of.interest = 2018,
                                                               legend.position.to.use = c(0.8, 0.8),
                                                               color.vector.to.use = color.vector.metabolism,
                                                               coverage_limits = c(0, 0.8),
                                                               depth_limits = c(40, 18),
                                                               taxonomy.column.name = "estimatedMetabolism",
                                                               ylab.to.use = "hgcA gene coverage")
pdf("results/hgcA_analysis/profiles/2018_profiles_286metalimnion_metabolism.pdf",
    height = 4,
    width = 4)
profile.286.2018.closeup
dev.off()

#### 2018 close-up profile of hypolimnion at 286, show inferred metabolism ####
profile.286.2018.closeup <- plot.profile.of.gene.with.taxonomy(data.to.use = hgcA.data,
                                                               RM.of.interest = 286,
                                                               year.of.interest = 2018,
                                                               legend.position.to.use = c(0.8, 0.8),
                                                               color.vector.to.use = color.vector.metabolism,
                                                               coverage_limits = c(0, 1),
                                                               depth_limits = c(80, 40),
                                                               taxonomy.column.name = "estimatedMetabolism",
                                                               ylab.to.use = "hgcA gene coverage")
pdf("results/hgcA_analysis/profiles/2018_profiles_286hypolimnion_metabolism.pdf",
    height = 4,
    width = 4)
profile.286.2018.closeup
dev.off()







#### 2019 profiles ####
profile.300.2019 <- plot.profile.of.gene.with.taxonomy(data.to.use = hgcA.data,
                                                       RM.of.interest = 300,
                                                       year.of.interest = 2019,
                                                       legend.position.to.use = c(0.6, 0.8),
                                                       depth_limits = c(60, 0),
                                                       coverage_limits = c(0, 0.8),
                                                       taxonomy.column.name = "manual_classification",
                                                       ylab.to.use = "hgcA gene coverage\n(per 100X SCG)",
                                                       titleToUse = "RM300 in 2019")
profile.310.2019 <- plot.profile.of.gene.with.taxonomy(data.to.use = hgcA.data,
                                                       RM.of.interest = 310,
                                                       year.of.interest = 2019,
                                                       legend.position.to.use = "none",
                                                       depth_limits = c(60, 0),
                                                       coverage_limits = c(0, 0.8),
                                                       taxonomy.column.name = "manual_classification",
                                                       ylab.to.use = "hgcA gene coverage\n(per 100X SCG)",
                                                       titleToUse = "RM310 in 2019")

pdf("results/hgcA_analysis/profiles/2019_profiles.pdf",
    height = 4.5,
    width = 6)
profile.300.2019 + profile.310.2019
dev.off()


#### 2019 close-up profile of hypolimnion at 300, show inferred metabolism ####
profile.300.2019 <- plot.profile.of.gene.with.taxonomy(data.to.use = hgcA.data,
                                                       RM.of.interest = 300,
                                                       year.of.interest = 2019,
                                                       legend.position.to.use = c(0.6, 0.8),
                                                       depth_limits = c(60, 0),
                                                       color.vector.to.use = color.vector.metabolism,
                                                       coverage_limits = c(0, 0.8),
                                                       taxonomy.column.name = "estimatedMetabolism",
                                                       ylab.to.use = "hgcA gene coverage\n(per 100X SCG)",
                                                       titleToUse = "RM300 in 2019")
profile.310.2019 <- plot.profile.of.gene.with.taxonomy(data.to.use = hgcA.data,
                                                       RM.of.interest = 310,
                                                       year.of.interest = 2019,
                                                       color.vector.to.use = color.vector.metabolism,
                                                       legend.position.to.use = "none",
                                                       depth_limits = c(60, 0),
                                                       coverage_limits = c(0, 0.8),
                                                       taxonomy.column.name = "estimatedMetabolism",
                                                       ylab.to.use = "hgcA gene coverage\n(per 100X SCG)",
                                                       titleToUse = "RM310 in 2019")
pdf("results/hgcA_analysis/profiles/2019_profiles_300hypolimnion_metabolism.pdf",
    height = 4,
    width = 4)
profile.300.2019 
dev.off()
pdf("results/hgcA_analysis/profiles/2019_profiles_310hypolimnion_metabolism.pdf",
    height = 4,
    width = 4)
profile.310.2019
dev.off()
