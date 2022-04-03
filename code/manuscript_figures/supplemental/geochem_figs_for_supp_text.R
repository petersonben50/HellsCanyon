#### code/manuscript_figures/geochem_figs_for_supp_text.R ####
# Benjamin D. Peterson


#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(naniar)
library(readxl)
library(patchwork)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/gene_plotting_functions.R")


#### Read in data ####
seabird.data <- readRDS("dataEdited/seabird/seabird_data.rds")
geochem.data <- read.csv("dataEdited/waterChemistry/geochem_WC.csv")
Hg.data <- read.csv("dataEdited/waterChemistry/Hg_data.csv") %>%
  group_by(RM, depth, date, constituent) %>%
  summarise(concentration = mean(concentration))


#### Sonde profile function ####
sonde.profile <- function(RM.of.interest,
                          year.of.interest,
                          depth.range.of.interest,
                          concentration.of.interest,
                          legend.location.of.interest) {
  seabird.data %>%
    mutate(DO = DO * 3) %>%
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
    scale_y_continuous(sec.axis = sec_axis(~ . / 3,
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
# Sonde profiles at 286 in 2019
sonde.286.2019 <- ggplot() + theme_void()

# Sonde profiles at 300 in 2019
sonde.300.2019 <- ggplot() + theme_void()

# Sonde profiles at 310 in 2019
sonde.310.2019 <- sonde.profile(310, 2019,
                                c(40, 0), c(0, 30),
                                "none")
# Sonde profile with legend
sonde.legend <- sonde.profile(310, 2019,
                              c(40, 0), c(0, 30),
                              c(0.25, 0.75))
# Check out plots
sonde.300.2019 + sonde.310.2019




#### Set redox plotting function ####
geochem.profile <- function(RM.of.interest,
                            year.of.interest,
                            depth.range.of.interest,
                            concentration.of.interest,
                            legend.location.of.interest) {
  geochem.data %>%
    filter(RM == RM.of.interest,
           year(date) == year.of.interest,
           month(date) == 7,
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
# Geochem profiles at 286 in 2019
geochem.286.2019 <- geochem.profile(286, 2019,
                                    c(80, 0), c(0, 2),
                                    "none")
# Geochem profiles at 300 in 2019
geochem.300.2019 <- geochem.profile(300, 2019,
                                    c(60, 0), c(0, 2),
                                    "none")
#c(0.7, 0.5)
# Geochem profiles at 310 in 2019
geochem.310.2019 <- geochem.profile(310, 2019,
                                    c(40, 0), c(0, 2),
                                    "none")
# Geochem profile with legend
geochem.legend <- geochem.profile(300, 2019,
                                  c(60, 0), c(0, 2),
                                  c(0.8, 0.8))

# Check out plots
(geochem.300.2019 + geochem.310.2019)




#### Set Hg plotting function ####
RM.of.interest <- 300
year.of.interest <- 2019
Hg.profile <- function(RM.of.interest,
                       year.of.interest,
                       depth.range.of.interest,
                       concentration.of.interest,
                       legend.location.of.interest) {
  Hg.data %>%
    filter(RM == RM.of.interest,
           year(date) == year.of.interest,
           month(date) == 7,
           day(date) >= 20) %>%
    spread(key = constituent,
           value = concentration) %>%
    mutate(FMHG_per = FMHG_per * (10/3)) %>%
    gather(key = constituent,
           value = concentration,
           -c(1:3)) %>%
    filter(!is.na(concentration),
           constituent %in% names(color.vector)) %>%
    mutate(date.site = paste(year(date), "-RM", RM,
                             sep = "")) %>%
    mutate(constituent = fct_relevel(constituent, names(color.vector))) %>%
    ggplot(aes(y = concentration,
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
# Hg profiles at 286 in 2019
Hg.2019.286 <- Hg.profile(286, 2019,
                          depth.range.of.interest = c(80, 0),
                          concentration.of.interest = c(0, 2),
                          legend.location.of.interest = "none")

# Hg profiles at 300 in 2019
Hg.2019.300 <- Hg.profile(300, 2019,
                          depth.range.of.interest = c(60, 0),
                          concentration.of.interest = c(0, 2),
                          legend.location.of.interest = "none")
# Hg profiles at 310 in 2019
Hg.2019.310 <- Hg.profile(310, 2019,
                          depth.range.of.interest = c(40, 0),
                          concentration.of.interest = c(0, 2),
                          legend.location.of.interest = "none")
# Hg profile with legend
Hg.legend <- Hg.profile(310, 2019,
                        depth.range.of.interest = c(40, 0),
                        concentration.of.interest = c(0, 2),
                        c(0.8, 0.8))
# Check out plots
Hg.2019.286 + Hg.2019.300 + Hg.2019.310




#### Save out needed plots ####
# RM 286 in 2019
pdf("results/manuscript_figures/geochem_figs_for_supp_text/RM286_2019.pdf",
    width = 5,
    height = 3.5)
sonde.286.2019 + geochem.286.2019 + Hg.2019.286
dev.off()
# RM 300 in 2019
pdf("results/manuscript_figures/geochem_figs_for_supp_text/RM300_2019.pdf",
    width = 5,
    height = 3.5)
sonde.300.2019 + geochem.300.2019 + Hg.2019.300
dev.off()
pdf("results/manuscript_figures/geochem_figs_for_supp_text/RM310_2019.pdf",
    width = 5,
    height = 3.5)
sonde.310.2019 + geochem.310.2019 + Hg.2019.310
dev.off()



#### Get fig with legends ####
pdf("results/manuscript_figures/geochem_figs_for_supp_text/for_legend_extraction.pdf",
    width = 5,
    height = 3.5)
sonde.legend + geochem.legend + Hg.legend
dev.off()
