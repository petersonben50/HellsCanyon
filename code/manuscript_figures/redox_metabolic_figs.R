#### code/manuscript_figures/redox_metabolic_figs.R ####
# Benjamin D. Peterson

#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(naniar)
library(readxl)
library(patchwork)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/gene_plotting_functions.R")


#### Read in geochem data ####
geochem.data.2017.2018 <- read.csv("dataEdited/waterChemistry/geochem_WC_2015_2018.csv") %>%
  filter(year(date) %in% c("2017", "2018")) %>%
  filter(month(date) == "9") %>%
  arrange(depth) %>%
  select(date, RM, depth, Fe.part_ug.L, Mn.part_ug.L,
         Mn_ug.L, NO3_mgN.L, sulfide_mg.L)
# Replace non-detects with zeroes in sulfide column
geochem.data.2017.2018$sulfide_mg.L[geochem.data.2017.2018$sulfide_mg.L == "nd"] <- 0
geochem.data.2017.2018$sulfide_mg.L <- as.numeric(geochem.data.2017.2018$sulfide_mg.L)

geochem.data.2019 <- read.csv("dataEdited/waterChemistry/geochem_WC_2019.csv") %>%
  select(date, RM, depth, Fe.part_ug.L, Mn.part_ug.L,
         Mn_ug.L, NO3_mgN.L)
geochem.data.2019$sulfide_mg.L <- NA

geochem.data <- rbind(geochem.data.2017.2018,
                      geochem.data.2019)
rm(geochem.data.2017.2018,
   geochem.data.2019)


#### Read in SeaBird data and combine with geochem ####
seabird.data <- readRDS("dataEdited/seabird/seabird_data.rds")


#### Read in metabolic gene data ####
teap.data <- readRDS("dataEdited/metabolic_analyses/depth/metabolicProtein_depth_clean.rds")
EET.data <- readRDS("dataEdited/metabolic_analyses/BBOMP/bbomp_depth_clean.rds")




#########################################
#### 2017 and 2018 profiles at RM286 ####


#### Sonde data ####
color.vector <- c(cb.translator["blue"],
                  cb.translator["black"])
names(color.vector) <- c("temp", "DO")
line.vector <- c(1, 3)
names(line.vector) <- c("temp", "DO")
labels.vector <- c("Temperature", "Dissolved O2")
names(labels.vector) <- c("temp", "DO")

# Sonde profiles at 286 in 2017
sonde.286.2017 <- seabird.data %>%
  mutate(DO = DO * 5) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:3)) %>%
  filter(RM == 286,
         year(date) == 2017,
         constituent %in% names(color.vector)) %>%
  ggplot(aes(x = depth,
             y = concentration,
             color = constituent)) +
  geom_point(aes(color = constituent)) +
  geom_line(aes(color = constituent,
                linetype = constituent)) +
  scale_colour_manual(values = color.vector,
                      labels = labels.vector) +
  scale_linetype_manual(values = line.vector,
                        labels = labels.vector) +
  coord_flip(xlim = c(75, 0),
             ylim = c(0, 30)) +
  ylab("Temp. (C)") +
  xlab("Depth (m)") +
  theme_classic() +
  theme(legend.position = c(0.75, 0.25),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) + 
  scale_y_continuous(sec.axis = sec_axis(~ . / 5,
                                         name = "DO (mg/L)"))
# Sonde profiles at 286 in 2018
sonde.286.2018 <- seabird.data %>%
  mutate(DO = DO * 5) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:3)) %>%
  filter(RM == 286,
         year(date) == 2018,
         constituent %in% names(color.vector)) %>%
  ggplot(aes(x = depth,
             y = concentration,
             color = constituent)) +
  geom_point(aes(color = constituent)) +
  geom_line(aes(color = constituent,
                linetype = constituent)) +
  scale_colour_manual(values = color.vector,
                      labels = labels.vector) +
  scale_linetype_manual(values = line.vector,
                        labels = labels.vector) +
  coord_flip(xlim = c(75, 0),
             ylim = c(0, 30)) +
  ylab("Temp. (C)") +
  xlab("Depth (m)") +
  theme_classic() +
  theme(legend.position = c(0.75, 0.25),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) + 
  scale_y_continuous(sec.axis = sec_axis(~ . / 5,
                                         name = "DO (mg/L)"))


sonde.286.2017 + sonde.286.2018

#### Low redox TEAs ####
color.vector <- c(cb.translator["skyblue"],
                  cb.translator["reddishpurple"],
                  cb.translator["vermillion"],
                  cb.translator["bluishgreen"])
names(color.vector) <- c("NO3_mgN.L", "Mn.part_mg.L",
                         "Mn_mg.L", "sulfide_mg.L")
line.vector <- c(1, 4, 3, 2)
names(line.vector) <- c("NO3_mgN.L", "Mn.part_mg.L",
                         "Mn_mg.L", "sulfide_mg.L")
labels.vector <- c("NO3 (mgN/L)", "Part. Mn",
                   "Diss. Mn", "Sulfide")
names(labels.vector) <- c("NO3_mgN.L", "Mn.part_mg.L",
                          "Mn_mg.L", "sulfide_mg.L")
points.vector <- c(16, 2, 17, 18)
names(points.vector) <- c("NO3_mgN.L", "Mn.part_mg.L",
                          "Mn_mg.L", "sulfide_mg.L")
# Low redox constituents at 286 in 2017
geochem.286.2017 <- geochem.data %>%
  mutate(Mn_mg.L = Mn_ug.L / 1000,
         Mn.part_mg.L = Mn.part_ug.L / 1000) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:3)) %>%
  filter(RM == 286,
         year(date) == 2017,
         constituent %in% names(color.vector)) %>%
  ggplot(aes(x = depth,
             y = concentration,
             color = constituent,
             linetype = constituent,
             shape = constituent)) +
  geom_point(aes(color = constituent)) +
  geom_line(aes(color = constituent)) +
  scale_colour_manual(values = color.vector,
                      labels = labels.vector) +
  scale_shape_manual(values = points.vector,
                     labels = labels.vector) +
  scale_linetype_manual(values = line.vector,
                        labels = labels.vector) +
  coord_flip(xlim = c(75, 0),
             ylim = c(0, 1.8)) +
  ylab("Concentration (mg/L)") +
  xlab("Depth (m)") +
  theme_classic() +
  theme(legend.position = c(0.75, 0.4),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
# Low redox constituents at 286 in 2018
geochem.286.2018 <- geochem.data %>%
  mutate(Mn_mg.L = Mn_ug.L / 1000,
         Mn.part_mg.L = Mn.part_ug.L / 1000) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:3)) %>%
  filter(RM == 286,
         year(date) == 2018,
         constituent %in% names(color.vector)) %>%
  ggplot(aes(x = depth,
             y = concentration,
             color = constituent,
             linetype = constituent,
             shape = constituent)) +
  geom_point(aes(color = constituent)) +
  geom_line(aes(color = constituent)) +
  scale_colour_manual(values = color.vector,
                      labels = labels.vector) +
  scale_shape_manual(values = points.vector,
                     labels = labels.vector) +
  scale_linetype_manual(values = line.vector,
                        labels = labels.vector) +
  coord_flip(xlim = c(75, 0),
             ylim = c(0, 1.8)) +
  ylab("Concentration (mg/L)") +
  xlab("Depth (m)") +
  theme_classic() +
  theme(legend.position = c(0.775, 0.8),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))



#### N cycling genes ####
color.vector <- c(cb.translator["vermillion"],
                  cb.translator["reddishpurple"],
                  cb.translator["bluishgreen"])
names(color.vector) <- c("narG", "napA", "nirS")
points.vector <- c(16, 2, 17)
names(points.vector) <- c("narG", "napA", "nirS")
N.cycling.286.2017 <- plot.profile.for.multiple.genes(marker.depth.df = teap.data,
                                                      genesOfInterest = names(color.vector),
                                                      yearOfInterest = "2017",
                                                      RMofInterest = "286",
                                                      show.mean.coverage = FALSE,
                                                      depth_limits = c(75, 0),
                                                      color.vector.to.use = color.vector,
                                                      point.vector.to.use = points.vector,
                                                      coverage_limits = c(0, 125),
                                                      titleToUse = element_blank(),
                                                      legend.position.to.use = c(0.7, 0.7))
N.cycling.286.2018 <- plot.profile.for.multiple.genes(marker.depth.df = teap.data,
                                                      genesOfInterest = names(color.vector),
                                                      yearOfInterest = "2018",
                                                      RMofInterest = "286",
                                                      show.mean.coverage = FALSE,
                                                      depth_limits = c(75, 0),
                                                      color.vector.to.use = color.vector,
                                                      point.vector.to.use = points.vector,
                                                      coverage_limits = c(0, 125),
                                                      titleToUse = element_blank(),
                                                      legend.position.to.use = c(0.7, 0.7))
N.cycling.286.2017 + N.cycling.286.2018


#### EET genes ####
color.vector <- c(cb.translator["orange"],
                  cb.translator["yellow"],
                  cb.translator["reddishpurple"])
names(color.vector) <- c("cluster_1_ExtE", "cluster_2_Omb", "cluster_3")
points.vector <- c(16, 2, 17)
names(points.vector) <- c("cluster_1_ExtE", "cluster_2_Omb", "cluster_3")
EET.286.2017 <- plot.profile.for.multiple.genes(marker.depth.df = EET.data,
                                                genesOfInterest = names(color.vector),
                                                yearOfInterest = "2017",
                                                RMofInterest = "286",
                                                gene.name.column = "classification",
                                                show.mean.coverage = FALSE,
                                                depth_limits = c(75, 0),
                                                color.vector.to.use = color.vector,
                                                point.vector.to.use = points.vector,
                                                coverage_limits = c(0, 1),
                                                titleToUse = element_blank(),
                                                legend.position.to.use = c(0.6, 0.5))
EET.286.2018 <- plot.profile.for.multiple.genes(marker.depth.df = EET.data,
                                                genesOfInterest = names(color.vector),
                                                yearOfInterest = "2018",
                                                RMofInterest = "286",
                                                gene.name.column = "classification",
                                                show.mean.coverage = FALSE,
                                                depth_limits = c(75, 0),
                                                color.vector.to.use = color.vector,
                                                point.vector.to.use = points.vector,
                                                coverage_limits = c(0, 1),
                                                titleToUse = element_blank(),
                                                legend.position.to.use = c(0.6, 0.6))
EET.286.2017 + EET.286.2018


#### S cycling genes ####
color.vector <- c(cb.translator["orange"],
                  cb.translator["reddishpurple"],
                  cb.translator["bluishgreen"],
                  cb.translator["skyblue"])
names(color.vector) <- c("rdsrA", "soxB", "dsrA", "dsrD")
points.vector <- c(16, 2, 17, 1)
names(points.vector) <- c("rdsrA", "soxB", "dsrA", "dsrD")
S.cycling.286.2017 <- plot.profile.for.multiple.genes(marker.depth.df = teap.data,
                                                      genesOfInterest = names(color.vector),
                                                      yearOfInterest = "2017",
                                                      RMofInterest = "286",
                                                      show.mean.coverage = FALSE,
                                                      depth_limits = c(75, 0),
                                                      color.vector.to.use = color.vector,
                                                      point.vector.to.use = points.vector,
                                                      coverage_limits = c(0, 90),
                                                      titleToUse = element_blank(),
                                                      legend.position.to.use = c(0.8, 0.8))
S.cycling.286.2018 <- plot.profile.for.multiple.genes(marker.depth.df = teap.data,
                                                      genesOfInterest = names(color.vector),
                                                      yearOfInterest = "2018",
                                                      RMofInterest = "286",
                                                      show.mean.coverage = FALSE,
                                                      depth_limits = c(75, 0),
                                                      color.vector.to.use = color.vector,
                                                      point.vector.to.use = points.vector,
                                                      coverage_limits = c(0, 90),
                                                      titleToUse = element_blank(),
                                                      legend.position.to.use = c(0.8, 0.8))
S.cycling.286.2017 + S.cycling.286.2018



#### Arrange plots ####
pdf("results/manuscript_figures/geochem_teap_RM286_2017_2018.pdf",
    width = 15,
    height = 9)
wrap_plots(sonde.286.2017, geochem.286.2017, N.cycling.286.2017, EET.286.2017, S.cycling.286.2017,
           sonde.286.2018, geochem.286.2018, N.cycling.286.2018, EET.286.2018, S.cycling.286.2018,
           nrow = 2)
dev.off()























#########################################
#### 2017 and 2018 profiles at RM300 ####


#### Sonde data ####
color.vector <- c(cb.translator["blue"],
                  cb.translator["black"])
names(color.vector) <- c("temp", "DO")
line.vector <- c(1, 3)
names(line.vector) <- c("temp", "DO")
labels.vector <- c("Temperature", "Dissolved O2")
names(labels.vector) <- c("temp", "DO")

# Sonde profiles at 300 in 2017
sonde.300.2017 <- seabird.data %>%
  mutate(DO = DO * 5) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:3)) %>%
  filter(RM == 300,
         year(date) == 2017,
         constituent %in% names(color.vector)) %>%
  ggplot(aes(x = depth,
             y = concentration,
             color = constituent)) +
  geom_point(aes(color = constituent)) +
  geom_line(aes(color = constituent,
                linetype = constituent)) +
  scale_colour_manual(values = color.vector,
                      labels = labels.vector) +
  scale_linetype_manual(values = line.vector,
                        labels = labels.vector) +
  coord_flip(xlim = c(60, 0),
             ylim = c(0, 35)) +
  ylab("Temp. (C)") +
  xlab("Depth (m)") +
  theme_classic() +
  theme(legend.position = c(0.75, 0.2),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) + 
  scale_y_continuous(sec.axis = sec_axis(~ . / 5,
                                         name = "DO (mg/L)"))
# Sonde profiles at 300 in 2018
sonde.300.2018 <- seabird.data %>%
  mutate(DO = DO * 5) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:3)) %>%
  filter(RM == 300,
         year(date) == 2018,
         constituent %in% names(color.vector)) %>%
  ggplot(aes(x = depth,
             y = concentration,
             color = constituent)) +
  geom_point(aes(color = constituent)) +
  geom_line(aes(color = constituent,
                linetype = constituent)) +
  scale_colour_manual(values = color.vector,
                      labels = labels.vector) +
  scale_linetype_manual(values = line.vector,
                        labels = labels.vector) +
  coord_flip(xlim = c(60, 0),
             ylim = c(0, 35)) +
  ylab("Temp. (C)") +
  xlab("Depth (m)") +
  theme_classic() +
  theme(legend.position = c(0.75, 0.2),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) + 
  scale_y_continuous(sec.axis = sec_axis(~ . / 5,
                                         name = "DO (mg/L)"))


sonde.300.2017 + sonde.300.2018

#### Low redox TEAs ####
color.vector <- c(cb.translator["skyblue"],
                  cb.translator["reddishpurple"],
                  cb.translator["vermillion"],
                  cb.translator["bluishgreen"])
names(color.vector) <- c("NO3_mgN.L", "Mn.part_mg.L",
                         "Mn_mg.L", "sulfide_mg.L")
line.vector <- c(1, 4, 3, 2)
names(line.vector) <- c("NO3_mgN.L", "Mn.part_mg.L",
                        "Mn_mg.L", "sulfide_mg.L")
labels.vector <- c("NO3 (mgN/L)", "Part. Mn",
                   "Diss. Mn", "Sulfide")
names(labels.vector) <- c("NO3_mgN.L", "Mn.part_mg.L",
                          "Mn_mg.L", "sulfide_mg.L")
points.vector <- c(16, 2, 17, 18)
names(points.vector) <- c("NO3_mgN.L", "Mn.part_mg.L",
                          "Mn_mg.L", "sulfide_mg.L")
# Low redox constituents at 300 in 2017
geochem.300.2017 <- geochem.data %>%
  mutate(Mn_mg.L = Mn_ug.L / 1000,
         Mn.part_mg.L = Mn.part_ug.L / 1000) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:3)) %>%
  filter(RM == 300,
         year(date) == 2017,
         day(date) > 20,
         constituent %in% names(color.vector)) %>%
  ggplot(aes(x = depth,
             y = concentration,
             color = constituent,
             linetype = constituent,
             shape = constituent)) +
  geom_point(aes(color = constituent)) +
  geom_line(aes(color = constituent)) +
  scale_colour_manual(values = color.vector,
                      labels = labels.vector) +
  scale_shape_manual(values = points.vector,
                     labels = labels.vector) +
  scale_linetype_manual(values = line.vector,
                        labels = labels.vector) +
  coord_flip(xlim = c(60, 0),
             ylim = c(0, 1.8)) +
  ylab("Concentration (mg/L)") +
  xlab("Depth (m)") +
  theme_classic() +
  theme(legend.position = c(0.3, 0.8),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
# Low redox constituents at 300 in 2018
geochem.300.2018 <- geochem.data %>%
  mutate(Mn_mg.L = Mn_ug.L / 1000,
         Mn.part_mg.L = Mn.part_ug.L / 1000) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:3)) %>%
  filter(RM == 300,
         year(date) == 2018,
         day(date) > 20,
         constituent %in% names(color.vector)) %>%
  ggplot(aes(x = depth,
             y = concentration,
             color = constituent,
             linetype = constituent,
             shape = constituent)) +
  geom_point(aes(color = constituent)) +
  geom_line(aes(color = constituent)) +
  scale_colour_manual(values = color.vector,
                      labels = labels.vector) +
  scale_shape_manual(values = points.vector,
                     labels = labels.vector) +
  scale_linetype_manual(values = line.vector,
                        labels = labels.vector) +
  coord_flip(xlim = c(60, 0),
             ylim = c(0, 1.8)) +
  ylab("Concentration (mg/L)") +
  xlab("Depth (m)") +
  theme_classic() +
  theme(legend.position = c(0.778, 0.8),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))



#### N cycling genes ####
color.vector <- c(cb.translator["vermillion"],
                  cb.translator["reddishpurple"],
                  cb.translator["bluishgreen"])
names(color.vector) <- c("narG", "napA", "nirS")
points.vector <- c(16, 2, 17)
names(points.vector) <- c("narG", "napA", "nirS")
N.cycling.300.2017 <- plot.profile.for.multiple.genes(marker.depth.df = teap.data,
                                                      genesOfInterest = names(color.vector),
                                                      yearOfInterest = "2017",
                                                      RMofInterest = "300",
                                                      show.mean.coverage = FALSE,
                                                      depth_limits = c(60, 0),
                                                      color.vector.to.use = color.vector,
                                                      point.vector.to.use = points.vector,
                                                      coverage_limits = c(0, 220),
                                                      titleToUse = element_blank(),
                                                      legend.position.to.use = c(0.7, 0.7))
N.cycling.300.2018 <- plot.profile.for.multiple.genes(marker.depth.df = teap.data,
                                                      genesOfInterest = names(color.vector),
                                                      yearOfInterest = "2018",
                                                      RMofInterest = "300",
                                                      show.mean.coverage = FALSE,
                                                      depth_limits = c(60, 0),
                                                      color.vector.to.use = color.vector,
                                                      point.vector.to.use = points.vector,
                                                      coverage_limits = c(0, 220),
                                                      titleToUse = element_blank(),
                                                      legend.position.to.use = c(0.7, 0.7))
N.cycling.300.2017 + N.cycling.300.2018


#### EET genes ####
color.vector <- c(cb.translator["orange"],
                  cb.translator["yellow"],
                  cb.translator["reddishpurple"])
names(color.vector) <- c("cluster_1_ExtE", "cluster_2_Omb", "cluster_3")
points.vector <- c(16, 2, 17)
names(points.vector) <- c("cluster_1_ExtE", "cluster_2_Omb", "cluster_3")
EET.300.2017 <- plot.profile.for.multiple.genes(marker.depth.df = EET.data,
                                                genesOfInterest = names(color.vector),
                                                yearOfInterest = "2017",
                                                RMofInterest = "300",
                                                gene.name.column = "classification",
                                                show.mean.coverage = FALSE,
                                                depth_limits = c(60, 0),
                                                color.vector.to.use = color.vector,
                                                point.vector.to.use = points.vector,
                                                coverage_limits = c(0, 1),
                                                titleToUse = element_blank(),
                                                legend.position.to.use = c(0.6, 0.5))
EET.300.2018 <- plot.profile.for.multiple.genes(marker.depth.df = EET.data,
                                                genesOfInterest = names(color.vector),
                                                yearOfInterest = "2018",
                                                RMofInterest = "300",
                                                gene.name.column = "classification",
                                                show.mean.coverage = FALSE,
                                                depth_limits = c(60, 0),
                                                color.vector.to.use = color.vector,
                                                point.vector.to.use = points.vector,
                                                coverage_limits = c(0, 1),
                                                titleToUse = element_blank(),
                                                legend.position.to.use = c(0.6, 0.6))
EET.300.2017 + EET.300.2018


#### S cycling genes ####
color.vector <- c(cb.translator["orange"],
                  cb.translator["bluishgreen"],
                  cb.translator["skyblue"])
names(color.vector) <- c("rdsrA", "dsrA", "dsrD")
points.vector <- c(16, 17, 1)
names(points.vector) <- c("rdsrA", "dsrA", "dsrD")
S.cycling.300.2017 <- plot.profile.for.multiple.genes(marker.depth.df = teap.data,
                                                      genesOfInterest = names(color.vector),
                                                      yearOfInterest = "2017",
                                                      RMofInterest = "300",
                                                      show.mean.coverage = FALSE,
                                                      depth_limits = c(60, 0),
                                                      color.vector.to.use = color.vector,
                                                      point.vector.to.use = points.vector,
                                                      coverage_limits = c(0, 50),
                                                      titleToUse = element_blank(),
                                                      legend.position.to.use = c(0.8, 0.8))
S.cycling.300.2018 <- plot.profile.for.multiple.genes(marker.depth.df = teap.data,
                                                      genesOfInterest = names(color.vector),
                                                      yearOfInterest = "2018",
                                                      RMofInterest = "300",
                                                      show.mean.coverage = FALSE,
                                                      depth_limits = c(60, 0),
                                                      color.vector.to.use = color.vector,
                                                      point.vector.to.use = points.vector,
                                                      coverage_limits = c(0, 50),
                                                      titleToUse = element_blank(),
                                                      legend.position.to.use = c(0.8, 0.8))
S.cycling.300.2017 + S.cycling.300.2018



#### Arrange plots ####
pdf("results/manuscript_figures/geochem_teap_RM300_2017_2018.pdf",
    width = 15,
    height = 10)
wrap_plots(sonde.300.2017, geochem.300.2017, N.cycling.300.2017, EET.300.2017, S.cycling.300.2017,
           sonde.300.2018, geochem.300.2018, N.cycling.300.2018, EET.300.2018, S.cycling.300.2018,
           nrow = 2)
dev.off()





















#########################################
#### 2019 profiles at RM300 and RM310 ####


#### Sonde data ####
color.vector <- c(cb.translator["blue"],
                  cb.translator["black"])
names(color.vector) <- c("temp", "DO")
line.vector <- c(1, 3)
names(line.vector) <- c("temp", "DO")
labels.vector <- c("Temperature", "Dissolved O2")
names(labels.vector) <- c("temp", "DO")

# Sonde profiles at 300 in 2017
sonde.300.2019 <- ggplot() + theme_void()
# Sonde profiles at 310 in 2019
sonde.310.2019 <- seabird.data %>%
  mutate(DO = DO * 2) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:3)) %>%
  filter(RM == 310,
         year(date) == 2019,
         constituent %in% names(color.vector)) %>%
  ggplot(aes(x = depth,
             y = concentration,
             color = constituent)) +
  geom_point(aes(color = constituent)) +
  geom_line(aes(color = constituent,
                linetype = constituent)) +
  scale_colour_manual(values = color.vector,
                      labels = labels.vector) +
  scale_linetype_manual(values = line.vector,
                        labels = labels.vector) +
  coord_flip(xlim = c(40, 0),
             ylim = c(0, 33)) +
  ylab("Temp. (C)") +
  xlab("Depth (m)") +
  theme_classic() +
  theme(legend.position = c(0.75, 0.2),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) + 
  scale_y_continuous(sec.axis = sec_axis(~ . / 2,
                                         name = "DO (mg/L)"))


sonde.300.2019 + sonde.310.2019

#### TEAs ####
color.vector <- c(cb.translator["skyblue"],
                  cb.translator["reddishpurple"],
                  cb.translator["vermillion"],
                  cb.translator["bluishgreen"])
names(color.vector) <- c("NO3_mgN.L", "Mn.part_mg.L",
                         "Mn_mg.L", "sulfide_mg.L")
line.vector <- c(1, 4, 3, 2)
names(line.vector) <- c("NO3_mgN.L", "Mn.part_mg.L",
                        "Mn_mg.L", "sulfide_mg.L")
labels.vector <- c("NO3 (mgN/L)", "Part. Mn",
                   "Diss. Mn", "Sulfide")
names(labels.vector) <- c("NO3_mgN.L", "Mn.part_mg.L",
                          "Mn_mg.L", "sulfide_mg.L")
points.vector <- c(16, 2, 17, 18)
names(points.vector) <- c("NO3_mgN.L", "Mn.part_mg.L",
                          "Mn_mg.L", "sulfide_mg.L")
# Redox constituents at 300 in 2019
geochem.300.2019 <- geochem.data %>%
  mutate(Mn_mg.L = Mn_ug.L / 1000,
         Mn.part_mg.L = Mn.part_ug.L / 1000) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:3)) %>%
  filter(RM == 300,
         year(date) == 2019,
         day(date) > 20,
         constituent %in% names(color.vector)) %>%
  ggplot(aes(x = depth,
             y = concentration,
             color = constituent,
             linetype = constituent,
             shape = constituent)) +
  geom_point(aes(color = constituent)) +
  geom_line(aes(color = constituent)) +
  scale_colour_manual(values = color.vector,
                      labels = labels.vector) +
  scale_shape_manual(values = points.vector,
                     labels = labels.vector) +
  scale_linetype_manual(values = line.vector,
                        labels = labels.vector) +
  coord_flip(xlim = c(60, 0),
             ylim = c(0, 1.8)) +
  ylab("Concentration (mg/L)") +
  xlab("Depth (m)") +
  theme_classic() +
  theme(legend.position = c(0.7, 0.6),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
# Low redox constituents at 310 in 2019
geochem.310.2019 <- geochem.data %>%
  mutate(Mn_mg.L = Mn_ug.L / 1000,
         Mn.part_mg.L = Mn.part_ug.L / 1000) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:3)) %>%
  filter(RM == 310,
         year(date) == 2019,
         day(date) > 20,
         constituent %in% names(color.vector)) %>%
  ggplot(aes(x = depth,
             y = concentration,
             color = constituent,
             linetype = constituent,
             shape = constituent)) +
  geom_point(aes(color = constituent)) +
  geom_line(aes(color = constituent)) +
  scale_colour_manual(values = color.vector,
                      labels = labels.vector) +
  scale_shape_manual(values = points.vector,
                     labels = labels.vector) +
  scale_linetype_manual(values = line.vector,
                        labels = labels.vector) +
  coord_flip(xlim = c(40, 0),
             ylim = c(0, 1.8)) +
  ylab("Concentration (mg/L)") +
  xlab("Depth (m)") +
  theme_classic() +
  theme(legend.position = c(0.778, 0.4),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
geochem.300.2019 + geochem.310.2019


#### N cycling genes ####
color.vector <- c(cb.translator["vermillion"],
                  cb.translator["reddishpurple"],
                  cb.translator["bluishgreen"])
names(color.vector) <- c("narG", "napA", "nirS")
points.vector <- c(16, 2, 17)
names(points.vector) <- c("narG", "napA", "nirS")
N.cycling.300.2019 <- plot.profile.for.multiple.genes(marker.depth.df = teap.data,
                                                      genesOfInterest = names(color.vector),
                                                      yearOfInterest = "2019",
                                                      RMofInterest = "300",
                                                      show.mean.coverage = FALSE,
                                                      depth_limits = c(60, 0),
                                                      color.vector.to.use = color.vector,
                                                      point.vector.to.use = points.vector,
                                                      coverage_limits = c(0, 220),
                                                      titleToUse = element_blank(),
                                                      legend.position.to.use = c(0.7, 0.7))
N.cycling.310.2019 <- plot.profile.for.multiple.genes(marker.depth.df = teap.data,
                                                      genesOfInterest = names(color.vector),
                                                      yearOfInterest = "2019",
                                                      RMofInterest = "310",
                                                      show.mean.coverage = FALSE,
                                                      depth_limits = c(40, 0),
                                                      color.vector.to.use = color.vector,
                                                      point.vector.to.use = points.vector,
                                                      coverage_limits = c(0, 220),
                                                      titleToUse = element_blank(),
                                                      legend.position.to.use = c(0.7, 0.7))
N.cycling.300.2019 + N.cycling.310.2019


#### EET genes ####
color.vector <- c(cb.translator["orange"],
                  cb.translator["yellow"],
                  cb.translator["reddishpurple"])
names(color.vector) <- c("cluster_1_ExtE", "cluster_2_Omb", "cluster_3")
points.vector <- c(16, 2, 17)
names(points.vector) <- c("cluster_1_ExtE", "cluster_2_Omb", "cluster_3")
EET.300.2019 <- plot.profile.for.multiple.genes(marker.depth.df = EET.data,
                                                genesOfInterest = names(color.vector),
                                                yearOfInterest = "2019",
                                                RMofInterest = "300",
                                                gene.name.column = "classification",
                                                show.mean.coverage = FALSE,
                                                depth_limits = c(60, 0),
                                                color.vector.to.use = color.vector,
                                                point.vector.to.use = points.vector,
                                                coverage_limits = c(0, 2),
                                                titleToUse = element_blank(),
                                                legend.position.to.use = c(0.7, 0.7))
EET.310.2019 <- plot.profile.for.multiple.genes(marker.depth.df = EET.data,
                                                genesOfInterest = names(color.vector),
                                                yearOfInterest = "2019",
                                                RMofInterest = "310",
                                                gene.name.column = "classification",
                                                show.mean.coverage = FALSE,
                                                depth_limits = c(40, 0),
                                                color.vector.to.use = color.vector,
                                                point.vector.to.use = points.vector,
                                                coverage_limits = c(0, 2),
                                                titleToUse = element_blank(),
                                                legend.position.to.use = c(0.7, 0.7))
EET.300.2019 + EET.310.2019


#### S cycling genes ####
color.vector <- c(cb.translator["orange"],
                  cb.translator["reddishpurple"],
                  cb.translator["bluishgreen"],
                  cb.translator["skyblue"])
names(color.vector) <- c("rdsrA", "soxB", "dsrA", "dsrD")
points.vector <- c(16, 2, 17, 1)
names(points.vector) <- c("rdsrA", "soxB", "dsrA", "dsrD")
S.cycling.300.2019 <- plot.profile.for.multiple.genes(marker.depth.df = teap.data,
                                                      genesOfInterest = names(color.vector),
                                                      yearOfInterest = "2019",
                                                      RMofInterest = "300",
                                                      show.mean.coverage = FALSE,
                                                      depth_limits = c(60, 0),
                                                      color.vector.to.use = color.vector,
                                                      point.vector.to.use = points.vector,
                                                      coverage_limits = c(0, 150),
                                                      titleToUse = element_blank(),
                                                      legend.position.to.use = c(0.7, 0.7))
S.cycling.310.2019 <- plot.profile.for.multiple.genes(marker.depth.df = teap.data,
                                                      genesOfInterest = names(color.vector),
                                                      yearOfInterest = "2019",
                                                      RMofInterest = "310",
                                                      show.mean.coverage = FALSE,
                                                      depth_limits = c(40, 0),
                                                      color.vector.to.use = color.vector,
                                                      point.vector.to.use = points.vector,
                                                      coverage_limits = c(0, 150),
                                                      titleToUse = element_blank(),
                                                      legend.position.to.use = c(0.7, 0.7))
S.cycling.300.2019 + S.cycling.310.2019



#### Arrange plots ####
pdf("results/manuscript_figures/geochem_teap_RM300_RM310_2019.pdf",
    width = 15,
    height = 10)
wrap_plots(sonde.300.2019, geochem.300.2019, N.cycling.300.2019, EET.300.2019, S.cycling.300.2019,
           sonde.310.2019, geochem.310.2019, N.cycling.310.2019, EET.310.2019, S.cycling.310.2019,
           nrow = 2)
dev.off()

