#### code/geochem/all_profiles.R ####
# Benjamin D. Peterson

#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(readxl)
library(tidyverse)
source("code/HCC_plotting_needs.R")

#### Read in geochem data ####
geochem.data <- read.csv("dataEdited/geochem/geochem_WC.csv") %>%
  filter(depth != "SW") %>%
  mutate(depth = as.numeric(depth))
DO.data <- readRDS("dataEdited/seabird/seabird_data.rds") %>%
  filter(!is.na(diss_oxy_mg_per_l))


Mn.plotting.factor = 5
Mn.color = "orange"
nitrate.plotting.factor = 5
nitrate.color = "reddishpurple"
MeHg.plotting.factor = 2.5
MeHg.color = "bluishgreen"
DO.plotting.factor = 4
DO.color = "black"
min.elevation = 545
max.elevation = 632
size.of.points = 1.5
line.width = 2.5



#### Function to plot redox ####
redox.plot <- function(geochem.data.to.use = geochem.data,
                       DO.data.to.use = DO.data,
                       date.range =  c("2016-09-07", "2016-09-08"),
                       DO.date = NULL,
                       RMs.to.use = NULL,
                       Mn.color = "orange",
                       use.nitrate.chloride.ratios.YES.or.NO = "NO") {

  geochem.data.date <- geochem.data.to.use %>%
    spread(key = constituent,
           value = concentration) %>%
    filter(date >= as.Date(date.range[1]) &
             date <= as.Date(date.range[2]))

  if (is.null(RMs.to.use)){
    RMs.to.use <- sort(unique(geochem.data.date$RM))
  }

  for (RM.of.interest in RMs.to.use) {

    geochem.data.site <- geochem.data.date %>%
      filter(RM == RM.of.interest) %>%
      arrange(depth)
    date.of.sampling.at.RM <- geochem.data.site %>%
      select(date) %>%
      unlist(use.names = FALSE) %>%
      unique()


    #### Set up plot ####
    plot(x = NULL,
         y = NULL,
         xlab = "",
         ylab = "Elevation (m)",
         xlim = c(0, 11),
         ylim = c(min.elevation, max.elevation),
         xaxt = "n")

    title(main = paste(date.of.sampling.at.RM, ": RM", RM.of.interest, sep = ""),
          cex.main = 0.8)
    #### Add DO data ####
    if (is.null(DO.date)) {
      DO.data.date <- DO.data.to.use %>%
        filter(date == date.of.sampling.at.RM &
                 RM == RM.of.interest)
    } else {
      DO.data.date <- DO.data.to.use %>%
        filter(date %in% DO.date,
               RM == RM.of.interest)
    }

    points(x = DO.data.date$diss_oxy_mg_per_l*0.5,
          y = DO.data.date$elevation_m,
          col = "grey30",
          pch = 18)
    lines(x = DO.data.date$diss_oxy_mg_per_l*0.5,
          y = DO.data.date$elevation_m,
          col = "grey30",
          lwd = 1.5)
    # Add axis for DO
    axis(1,
         line = 6,
         at = seq(0, 10, by = 2.5),
         labels = seq(0, 10, by = 2.5)/0.5,
         cex.axis = 0.85)
    # Add title for DO
    title(xlab = "Dissolved O2 (mg/L)",
          line = 7.2,
          cex.lab = 0.9)
    
    #### Add Mn ####
    points(x = geochem.data.site$f_mn_mg_per_l*Mn.plotting.factor,
           y = geochem.data.site$elevation_m,
           col = cb.translator[Mn.color],
           pch = 18,
           cex = size.of.points)
    lines(x = geochem.data.site$f_mn_mg_per_l*Mn.plotting.factor,
          y = geochem.data.site$elevation_m,
          col = cb.translator[Mn.color],
          lwd = line.width)
    # Add axis for Mn and nitrate
    axis(1,
         at = seq(0, 10, by = 2.5),
         labels = seq(0, 10, by = 2.5)/Mn.plotting.factor,
         cex.axis = 0.85)
    # Add title for Mn and nitrate
    title(xlab = "Nitrate (mgN/L), Mn (mg/L)",
          line = 1.2,
          cex.lab = 0.9)
    
    
    #### Add nitrate ####
    points(x = geochem.data.site$f_no3_mg_n_per_l*nitrate.plotting.factor,
           y = geochem.data.site$elevation_m,
           col = cb.translator[nitrate.color],
           pch = 18,
           cex = size.of.points)
    lines(x = geochem.data.site$f_no3_mg_n_per_l*nitrate.plotting.factor,
          y = geochem.data.site$elevation_m,
          col = cb.translator[nitrate.color],
          lwd = line.width)
    
    #### Add MeHg ####
    points(x = geochem.data.site$MeHg_diss_ngL*MeHg.plotting.factor,
           y = geochem.data.site$elevation_m,
           col = cb.translator[MeHg.color],
           pch = 18,
           cex = size.of.points)
    lines(x = geochem.data.site$MeHg_diss_ngL*MeHg.plotting.factor,
          y = geochem.data.site$elevation_m,
          col = cb.translator[MeHg.color],
          lwd = line.width)
    # Add axis for MeHg
    axis(1,
         line = 3,
         at = seq(0, 10, by = 2.5),
         labels = seq(0, 10, by = 2.5)/MeHg.plotting.factor,
         cex.axis = 0.85)
    # Add title for MeHg
    title(xlab = "MeHg (ng/L)",
          line = 4.2,
          cex.lab = 0.9)
    
    
    }
  }



#### Save out plots ####
pdf("results/geochem/profiles_for_figure/WC_profiles_late_stratification.pdf",
    width = 7.5,
    height = 7)

par(mfrow = c(2,6),
    mar = c(9, 2.5, 1, 0.5),
    mgp=c(1.5,0.2,0),
    tck=-0.008)
redox.plot(date.range =  c("2016-10-03","2016-10-06"),
           RMs.to.use = c(286, 300, 310),
           DO.date = "2016-09-21")
redox.plot(date.range =  c("2017-09-25", "2017-09-28"),
           RMs.to.use = c(286, 300, 310))
redox.plot(date.range =  c("2018-09-24", "2018-09-26"),
           RMs.to.use = c(286, 300, 310))

dev.off()



#### Save out plots ####
pdf("results/geochem/profiles_for_figure/WC_profiles_early_stratification.pdf",
    width = 7.5,
    height = 7)

par(mfrow = c(2,6),
    mar = c(9, 2.5, 1, 0.5),
    mgp=c(1.5,0.2,0),
    tck=-0.008)
redox.plot(date.range = c("2017-06-05","2017-06-08"),
           RMs.to.use = c(286, 300, 310),
           DO.date = "2017-06-14")
redox.plot(date.range =  c("2018-06-18", "2018-06-19"),
           RMs.to.use = c(286, 300, 310),
           DO.date = "2018-06-12")
redox.plot(date.range =  c("2019-07-22", "2019-07-25"),
           RMs.to.use = c(286, 300, 310))
dev.off()

