#### code/manuscript_figures/geochem_figs_for_main_text.R ####
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


#### Read in data ####
seabird.data <- readRDS("dataEdited/seabird/seabird_data.rds")
geochem.data <- read.csv("dataEdited/waterChemistry/geochem_WC.csv")
Hg.data <- read.csv("dataEdited/waterChemistry/Hg_data.csv")


#### Sonde profile function ####
sonde.profile <- function(RM.of.interest,
                          year.of.interest,
                          depth.range.of.interest,
                          concentration.of.interest,
                          legend.location.of.interest) {
  seabird.data %>%
    mutate(DO = DO * 5) %>%
    gather(key = constituent,
           value = concentration,
           -c(1:3)) %>%
    filter(RM == RM.of.interest,
           year(date) == year.of.interest,
           constituent %in% names(color.vector)) %>%
    ggplot(aes(x = depth,
               y = concentration,
               color = constituent)) +
    geom_point(aes(color = constituent),
               size = 1) +
    geom_line(aes(color = constituent,
                  linetype = constituent)) +
    scale_colour_manual(values = color.vector,
                        labels = labels.vector) +
    scale_linetype_manual(values = line.vector,
                          labels = labels.vector) +
    coord_flip(xlim = depth.range.of.interest,
               ylim = concentration.of.interest) +
    ylab("Temp. (C)") +
    xlab("Depth (m)") +
    theme_classic() +
    theme(legend.position = legend.location.of.interest,
          legend.title = element_blank(),
          legend.key = element_rect(fill = "transparent", colour = "black"),
          legend.key.size = unit(1.75, 'lines'),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black")) + 
    scale_y_continuous(sec.axis = sec_axis(~ . / 5,
                                           name = "DO (mg/L)"))
}


#### Vectors for sonde profiles ####
color.vector <- c(cb.translator["skyblue"],
                  cb.translator["black"])
names(color.vector) <- c("temp", "DO")
line.vector <- c(1, 3)
names(line.vector) <- c("temp", "DO")
labels.vector <- c("Temperature", "Dissolved O2")
names(labels.vector) <- c("temp", "DO")


#### Generate sonde profiles ####
# Sonde profiles at 286 in 2017
sonde.2017.286 <- sonde.profile(286, 2017,
                                c(75, 0), c(0, 30),
                                "none")
# Sonde profiles at 286 in 2018
sonde.2018.286 <- sonde.profile(286, 2018,
                                c(75, 0), c(0, 35),
                                "none")
# Sonde profiles at 300 in 2017
sonde.2017.300 <- sonde.profile(300, 2017,
                                c(60, 0), c(0, 35),
                                "none")
# Sonde profiles at 300 in 2018
sonde.2018.300 <- sonde.profile(300, 2018,
                                c(60, 0), c(0, 35),
                                "none")
# Sonde profile with legend
sonde.legend <- sonde.profile(300, 2018,
                              c(60, 0), c(0, 35),
                              c(0.25, 0.75))
# Check out plots
sonde.2017.286 + sonde.2018.286 + sonde.2017.300 + sonde.2018.300




#### Set redox plotting function ####
geochem.profile <- function(RM.of.interest,
                            year.of.interest,
                            depth.range.of.interest,
                            concentration.of.interest,
                            legend.location.of.interest) {
  geochem.data %>%
    filter(RM == RM.of.interest,
           year(date) == year.of.interest,
           month(date) == 9,
           day(date) >= 20) %>%
    group_by(RM, depth, date, constituent) %>%
    summarize(concentration = mean(concentration)) %>%
    spread(key = constituent,
           value = concentration) %>%
    mutate(Mn_mg.L = Mn_ug.L / 1000,
           Mn.part_mg.L = Mn.part_ug.L / 1000) %>%
    gather(key = constituent,
           value = concentration,
           -c(1:3)) %>%
    filter(constituent %in% names(color.vector)) %>%
    ggplot(aes(x = depth,
               y = concentration,
               color = constituent,
               linetype = constituent,
               shape = constituent)) +
    geom_point(aes(color = constituent),
               size = 2) +
    geom_line(aes(color = constituent)) +
    scale_colour_manual(values = color.vector,
                        labels = labels.vector) +
    scale_shape_manual(values = points.vector,
                       labels = labels.vector) +
    scale_linetype_manual(values = line.vector,
                          labels = labels.vector) +
    coord_flip(xlim = depth.range.of.interest,
               ylim = concentration.of.interest) +
    ylab("Concentration (mg/L)") +
    xlab("Depth (m)") +
    theme_classic() +
    theme(legend.position = legend.location.of.interest,
          legend.title = element_blank(),
          legend.key = element_rect(fill = "transparent", colour = "black"),
          legend.key.size = unit(1.75, 'lines'),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"))
}


#### Vectors for geochem ####
color.vector <- c(cb.translator["reddishpurple"],
                  cb.translator["orange"],
                  cb.translator["blue"])
names(color.vector) <- c("NO3_mgN.L", "Mn_mg.L", "sulfide_mg.L")
line.vector <- c(1,3, 2)
names(line.vector) <- names(color.vector)
labels.vector <- c("NO3 (mgN/L)", "Diss. Mn", "Sulfide")
names(labels.vector) <- names(color.vector)
points.vector <- c(16, 17, 18)
names(points.vector) <- c("NO3_mgN.L", "Mn_mg.L", "sulfide_mg.L")


#### Generate geochem profiles ####
# Geochem profiles at 286 in 2017
geochem.2017.286 <- geochem.profile(286, 2017,
                                    c(75, 0), c(0, 2),
                                    "none")
#c(0.7, 0.5)
# Geochem profiles at 286 in 2018
geochem.2018.286 <- geochem.profile(286, 2018,
                                    c(75, 0), c(0, 2),
                                    "none")
#c(0.75, 0.75)
# Geochem profiles at 300 in 2017
geochem.2017.300 <- geochem.profile(300, 2017,
                                    c(60, 0), c(0, 2),
                                    "none")
#c(0.30, 0.75)
# Geochem profiles at 300 in 2018
geochem.2018.300 <- geochem.profile(300, 2018,
                                    c(60, 0), c(0, 2),
                                    "none")
#c(0.75, 0.75)

# Geochem profile with legend
geochem.legend <- geochem.profile(300, 2018,
                                  c(60, 0), c(0, 2),
                                  c(0.8, 0.8))

# Check out plots
(geochem.2017.286 + geochem.2018.286) / (geochem.2017.300 + geochem.2018.300)




#### Set Hg plotting function ####
RM.of.interest <- 286
year.of.interest <- 2017
Hg.profile <- function(RM.of.interest,
                       year.of.interest,
                       depth.range.of.interest,
                       concentration.of.interest,
                       legend.location.of.interest) {
  Hg.data %>%
    filter(RM == RM.of.interest,
           year(date) == year.of.interest,
           month(date) == 9,
           day(date) >= 20) %>%
    mutate(FMHG_per = FMHG_per * (10/3)) %>%
    gather(key = constituent,
           value = abundance,
           -c(1:3)) %>%
    filter(!is.na(abundance),
           constituent %in% names(color.vector)) %>%
    mutate(date.site = paste(year(date), "-RM", RM,
                             sep = "")) %>%
    mutate(constituent = fct_relevel(constituent, names(color.vector))) %>%
    ggplot(aes(y = abundance,
               x = depth,
               group = constituent,
               shape = constituent,
               linetype = constituent)) +
    geom_point(aes(color = constituent),
               size = 2) +
    geom_line(aes(color = constituent)) +
    scale_colour_manual(values = color.vector,
                        labels = labels.vector) +
    scale_shape_manual(values = points.vector,
                       labels = labels.vector) +
    scale_linetype_manual(values = line.vector,
                          labels = labels.vector) +
    coord_flip(xlim = depth.range.of.interest,
               ylim = concentration.of.interest) +
    theme_classic() +
    theme(legend.position = legend.location.of.interest,
          legend.title = element_blank(),
          legend.key = element_rect(fill = "transparent", colour = "black"),
          legend.key.size = unit(1.75, 'lines'),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black")) +
    xlab("Depth (m)") +
    ylab("Concentration (ng/L)") + 
    scale_y_continuous(sec.axis = sec_axis(~ . / (10/3),
                                           name = "Percent MeHg"))
}


#### Vectors for Hg ####
color.vector <- c(cb.translator["vermillion"],
                  cb.translator["bluishgreen"],
                  cb.translator["black"])
names(color.vector) <- c("FMHG", "FiHg", "FMHG_per")
line.vector <- c(1, 2, 4)
names(line.vector) <- names(color.vector)
labels.vector <- c("Filter-passing\nMeHg (ng/L)",
                   "Filter-passing\niHg (ng/L)",
                   "Filter-passing\nMeHg (%)")
names(labels.vector) <- names(color.vector)
points.vector <- c(2, 16, 17)
names(points.vector) <- names(color.vector)



#### Generate Hg profiles ####
# Hg profiles at 286 in 2017
Hg.2017.286 <- Hg.profile(286, 2017,
                          c(75, 0), c(0, 3.3),
                          "none")
# Hg profiles at 286 in 2018
Hg.2018.286 <- Hg.profile(286, 2018,
                          c(75, 0), c(0, 3.3),
                          "none")
# Hg profiles at 300 in 2017
Hg.2017.300 <- Hg.profile(300, 2017,
                          c(60, 0), c(0, 3.3),
                          "none")
# Hg profiles at 300 in 2018
Hg.2018.300 <- Hg.profile(300, 2018,
                          c(60, 0), c(0, 3.3),
                          "none")
# Hg profile with legend
Hg.legend <- Hg.profile(300, 2018,
                        c(60, 0), c(0, 3.3),
                        c(0.8, 0.8))
# Check out plots
(Hg.2017.286 + Hg.2018.286) / (Hg.2017.300 + Hg.2018.300)




#### Save out needed plots ####
# RM 286 in 2017
pdf("results/manuscript_figures/geochem_figs_for_main_text/RM286_2017.pdf",
    width = 7,
    height = 4)
sonde.2017.286 + geochem.2017.286 + Hg.2017.286
dev.off()
pdf("results/manuscript_figures/geochem_figs_for_main_text/RM300_2017.pdf",
    width = 7,
    height = 4)
sonde.2017.300 + geochem.2017.300 + Hg.2017.300
dev.off()

pdf("results/manuscript_figures/geochem_figs_for_main_text/RM286_2018.pdf",
    width = 7,
    height = 4)
sonde.2018.286 + geochem.2018.286 + Hg.2018.286
dev.off()

pdf("results/manuscript_figures/geochem_figs_for_main_text/RM300_2018.pdf",
    width = 7,
    height = 4)
sonde.2018.300 + geochem.2018.300 + Hg.2018.300
dev.off()



#### Get fig with legends ####

pdf("results/manuscript_figures/geochem_figs_for_main_text/for_legend_extraction.pdf",
    width = 7,
    height = 4)
sonde.legend + geochem.legend + Hg.legend
dev.off()
