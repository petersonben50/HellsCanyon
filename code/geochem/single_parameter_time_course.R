#### code/geochem/single_parameter_time_course.R ####
# Benjamin D. Peterson

#### Set the table ####

rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(ggpubr)
library(lubridate)
library(readxl)
library(tidyverse)
source("code/geochem/profile_functions.R")


#### Read in geochem data ####
geochem.data <- read.csv("dataEdited/geochem/geochem_WC.csv") %>%
  filter(depth != "SW") %>%
  mutate(depth = as.numeric(depth))
seabird.data <- readRDS("dataEdited/seabird/seabird_data.rds") %>%
  filter(!is.na(diss_oxy_mg_per_l))



#### Generate nitrate plots ####
nitrate.2015 <- time.course.profile.plot(geochem.data.to.use = geochem.data,
                         years.to.use = 2015,
                         RMs.to.use = c(286,300, 310),
                         parameter.to.plot = "f_no3_mg_n_per_l",
                         xlabel.to.use = "Nitrate (mgN/L)",
                         color.ramp.to.use = cb.translator[c("black", "orange", "yellow")]) +
    facet_wrap(~RM)
nitrate.2016 <- time.course.profile.plot(geochem.data.to.use = geochem.data,
                                         years.to.use = 2016,
                                         RMs.to.use = c(286,300, 310),
                                         parameter.to.plot = "f_no3_mg_n_per_l",
                                         xlabel.to.use = "Nitrate (mgN/L)",
                                         color.ramp.to.use = cb.translator[c("black", "orange", "yellow")]) +
  facet_wrap(~RM)

nitrate.2017 <- time.course.profile.plot(geochem.data.to.use = geochem.data,
                                         years.to.use = 2017,
                                         RMs.to.use = c(286,300, 310),
                                         parameter.to.plot = "f_no3_mg_n_per_l",
                                         xlabel.to.use = "Nitrate (mgN/L)",
                                         color.ramp.to.use = cb.translator[c("black", "orange", "yellow")]) +
  facet_wrap(~RM)
nitrate.2018 <- time.course.profile.plot(geochem.data.to.use = geochem.data,
                                         years.to.use = 2018,
                                         RMs.to.use = c(286,300, 310),
                                         parameter.to.plot = "f_no3_mg_n_per_l",
                                         xlabel.to.use = "Nitrate (mgN/L)",
                                         color.ramp.to.use = cb.translator[c("black", "orange", "yellow")]) +
  facet_wrap(~RM)
nitrate.2019 <- time.course.profile.plot(geochem.data.to.use = geochem.data,
                                         years.to.use = 2019,
                                         RMs.to.use = c(286,300, 310),
                                         parameter.to.plot = "f_no3_mg_n_per_l",
                                         xlabel.to.use = "Nitrate (mgN/L)",
                                         color.ramp.to.use = cb.translator[c("black", "orange", "yellow")]) +
  facet_wrap(~RM)


#### Save out nitrate plots ####
pdf("results/geochem/nitrate_time_course.pdf",
    width = 7,
    height = 20)
ggarrange(nitrate.2015, nitrate.2016, nitrate.2017, nitrate.2018, nitrate.2019,
          ncol = 1)
dev.off()
rm(nitrate.2015, nitrate.2016, nitrate.2017, nitrate.2018, nitrate.2019)



#### Generate MeHg plots ####
MeHg.2015 <- time.course.profile.plot(geochem.data.to.use = geochem.data,
                                      years.to.use = 2015,
                                      RMs.to.use = c(286,300, 310),
                                      parameter.to.plot = "MeHg_diss_ngL",
                                      xlabel.to.use = "MeHg (ng/L)",
                                      color.ramp.to.use = cb.translator[c("black", "orange", "yellow")],
                                      concentrations.to.use = c(0,4)) +
  facet_wrap(~RM)
MeHg.2016 <- time.course.profile.plot(geochem.data.to.use = geochem.data,
                                         years.to.use = 2016,
                                         RMs.to.use = c(286,300, 310),
                                         parameter.to.plot = "MeHg_diss_ngL",
                                         xlabel.to.use = "MeHg (ng/L)",
                                         color.ramp.to.use = cb.translator[c("black", "orange", "yellow")]) +
  facet_wrap(~RM)

MeHg.2017 <- time.course.profile.plot(geochem.data.to.use = geochem.data,
                                      years.to.use = 2017,
                                      RMs.to.use = c(286,300, 310),
                                      parameter.to.plot = "MeHg_diss_ngL",
                                      xlabel.to.use = "MeHg (ng/L)",
                                      color.ramp.to.use = cb.translator[c("black", "orange", "yellow")],
                                      concentrations.to.use = c(0,4)) +
  facet_wrap(~RM)
MeHg.2018 <- time.course.profile.plot(geochem.data.to.use = geochem.data,
                                         years.to.use = 2018,
                                         RMs.to.use = c(286,300, 310),
                                         parameter.to.plot = "MeHg_diss_ngL",
                                         xlabel.to.use = "MeHg (ng/L)",
                                         color.ramp.to.use = cb.translator[c("black", "orange", "yellow")]) +
  facet_wrap(~RM)
MeHg.2019 <- time.course.profile.plot(geochem.data.to.use = geochem.data,
                                         years.to.use = 2019,
                                         RMs.to.use = c(286,300, 310),
                                         parameter.to.plot = "MeHg_diss_ngL",
                                         xlabel.to.use = "MeHg (ng/L)",
                                         color.ramp.to.use = cb.translator[c("black", "orange", "yellow")]) +
  facet_wrap(~RM)


#### Save out MeHg plots ####
pdf("results/geochem/MeHg_time_course.pdf",
    width = 7,
    height = 20)
ggarrange(MeHg.2015, MeHg.2016, MeHg.2017, MeHg.2018, MeHg.2019,
          ncol = 1)
dev.off()
rm(MeHg.2015, MeHg.2016, MeHg.2017, MeHg.2018, MeHg.2019)






#### Generate DO plots ####
DO.2015 <- seabird.time.course.profile.plot(seabird.data.to.use = seabird.data,
                                            years.to.use = 2015,
                                            RMs.to.use = c(286,300, 310),
                                            parameter.to.plot = "diss_oxy_mg_per_l",
                                            xlabel.to.use = "DO (mg/L)",
                                            color.ramp.to.use = cb.translator[c("skyblue", "bluishgreen", "black")],
                                            concentrations.to.use = c(0,15)) +
  facet_wrap(~RM)
DO.2016 <- seabird.time.course.profile.plot(seabird.data.to.use = seabird.data,
                                            years.to.use = 2016,
                                            RMs.to.use = c(286,300, 310),
                                            parameter.to.plot = "diss_oxy_mg_per_l",
                                            xlabel.to.use = "DO (mg/L)",
                                            color.ramp.to.use = cb.translator[c("skyblue", "bluishgreen", "black")],
                                            concentrations.to.use = c(0,15)) +
  facet_wrap(~RM)

DO.2017 <- seabird.time.course.profile.plot(seabird.data.to.use = seabird.data,
                                            years.to.use = 2017,
                                            RMs.to.use = c(286,300, 310),
                                            parameter.to.plot = "diss_oxy_mg_per_l",
                                            xlabel.to.use = "DO (mg/L)",
                                            color.ramp.to.use = cb.translator[c("skyblue", "bluishgreen", "black")],
                                            concentrations.to.use = c(0,15)) +
  facet_wrap(~RM)
DO.2018 <- seabird.time.course.profile.plot(seabird.data.to.use = seabird.data,
                                            years.to.use = 2018,
                                            RMs.to.use = c(286,300, 310),
                                            parameter.to.plot = "diss_oxy_mg_per_l",
                                            xlabel.to.use = "DO (mg/L)",
                                            color.ramp.to.use = cb.translator[c("skyblue", "bluishgreen", "black")],
                                            concentrations.to.use = c(0,15)) +
  facet_wrap(~RM)
DO.2019 <- seabird.time.course.profile.plot(seabird.data.to.use = seabird.data,
                                            years.to.use = 2019,
                                            RMs.to.use = c(286,300, 310),
                                            parameter.to.plot = "diss_oxy_mg_per_l",
                                            xlabel.to.use = "DO (mg/L)",
                                            color.ramp.to.use = cb.translator[c("skyblue", "bluishgreen", "black")],
                                            concentrations.to.use = c(0,15)) +
  facet_wrap(~RM)


#### Save out DO plots ####
pdf("results/geochem/DO_time_course.pdf",
    width = 7,
    height = 20)
ggarrange(DO.2015, DO.2016, DO.2017, DO.2018, DO.2019,
          ncol = 1)
dev.off()
rm(DO.2015, DO.2016, DO.2017, DO.2018, DO.2019)
