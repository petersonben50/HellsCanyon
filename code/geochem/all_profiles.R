#### code/geochem/all_profiles.R ####
# Benjamin D. Peterson

#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(readxl)
library(tidyverse)
source("code/geochem/profile_functions.R")


#### Read in geochem data ####
geochem.data <- read.csv("dataEdited/waterChemistry/geochem_WC.csv") %>%
  filter(depth != "SW") %>%
  mutate(depth = as.numeric(depth))




#### List of sites ####
all.dates <- sort(unique(c(geochem.data$date)))


#### Plot temp/DO, redox, and MeHg in separate plots data ####

nitrate.plotting.factor.to.use = 5
sulfide.plotting.factor.to.use = 6+2/3
Mn.plotting.factor.to.use = 5
MeHg.plotting.factor.to.use = 2.5



# pdf("~/Downloads/test_geochem_plots.pdf",
#     width = 8.5,
#     height = 11)

#### 2014 profiles (won't use these, just HgT and MeHg) ####
all.dates[year(all.dates) == 2014]
# date.list <- list(
  # c("2014-04-23", "2014-04-23"),
  # c("2014-05-20", "2014-05-20"),
  # c("2014-06-17", "2014-06-17"),
  # c("2014-07-16", "2014-07-16"),
  # c("2014-08-12", "2014-08-13"),
  # c("2014-09-23", "2014-09-24"),
  # c("2014-10-14", "2014-10-15")
# )
# par(mfrow = c(4,4),
#     mar = c(9, 3, 1, 1),
#     mgp=c(1.5,0.4,0),
#     tck=-0.008)
# lapply(date.list,
#        function(list.vector) {
#          redox.plot(geochem.data.to.use = geochem.data,
#                     date.range = list.vector,
#                     RMs.to.use = c(286, 300, 310, 318),
#                     nitrate.plotting.factor = nitrate.plotting.factor.to.use,
#                     MeHg.plotting.factor = 2.5,
#                     plot.Mn.instead.of.sulfide.YES.or.NO = "YES",
#                     Mn.plotting.factor = Mn.plotting.factor.to.use,
#                     max.depth = 80)
#          
#        })




#### 2015 profiles: No sulfide data, plot Mn instead ####
all.dates[year(all.dates) == 2015]
date.list <- list(
  # c("2015-04-20", "2015-04-21"),
  # c("2015-05-18", "2015-05-19"),
  # c("2015-06-04", "2015-06-08"),
  c("2015-07-14", "2015-07-15"),
  c("2015-08-10", "2015-08-11"),
  c("2015-09-08", "2015-09-09")
  )
par(mfrow = c(3,4),
    mar = c(9, 3, 1, 1),
    mgp=c(1.5,0.4,0),
    tck=-0.008)
lapply(date.list,
       function(list.vector) {
         redox.plot(geochem.data.to.use = geochem.data,
                    date.range = list.vector,
                    RMs.to.use = c(286, 300, 310, 318),
                    nitrate.plotting.factor = nitrate.plotting.factor.to.use,
                    MeHg.plotting.factor = 2.5,
                    plot.Mn.instead.of.sulfide.YES.or.NO = "YES",
                    Mn.plotting.factor = Mn.plotting.factor.to.use,
                    max.depth = 80)
         
       })




#### 2016 profiles: Sulfide only measured in October, plot Mn for other dates instead ####
all.dates[year(all.dates) == 2016]
par(mfrow = c(4,4),
    mar = c(9, 3, 1, 1),
    mgp=c(1.5,0.4,0),
    tck=-0.008)

date.list <- list(
  # c("2016-04-05", "2016-04-06"),
  c("2016-05-17", "2016-05-18"),
  c("2016-07-12", "2016-07-13"),
  # c("2016-08-09", "2016-08-09"),
  c("2016-09-07", "2016-09-08")
  # c("2016-10-03","2016-10-06"),
  # c("2016-11-15", "2016-11-22")
  )
lapply(date.list,
       function(list.vector) {
         redox.plot(geochem.data.to.use = geochem.data,
                    date.range = list.vector,
                    RMs.to.use = c(286, 300, 310, 318),
                    nitrate.plotting.factor = nitrate.plotting.factor.to.use,
                    MeHg.plotting.factor = 2.5,
                    plot.Mn.instead.of.sulfide.YES.or.NO = "YES",
                    Mn.plotting.factor = Mn.plotting.factor.to.use,
                    max.depth = 80)
         
       })
redox.plot(geochem.data.to.use = geochem.data,
           date.range = c("2016-10-03","2016-10-06"),
           RMs.to.use = c(286, 300, 310, 318),
           nitrate.plotting.factor = nitrate.plotting.factor.to.use,
           MeHg.plotting.factor = 2.5,
           plot.Mn.instead.of.sulfide.YES.or.NO = "NO",
           sulfide.plotting.factor = sulfide.plotting.factor.to.use,
           max.depth = 80)



#### 2017 profiles:  ####
all.dates[year(all.dates) == 2017]
date.list <- list(
  # c("2017-03-20", "2017-03-21"),
  # c("2017-05-01", "2017-05-02"),
  c("2017-06-05", "2017-06-08"),
  c("2017-07-25", "2017-07-26"),
  c("2017-08-23", "2017-08-24"),
  c("2017-09-25", "2017-09-28")
  # c("2017-11-14", "2017-11-15")
  )
par(mfrow = c(4,4),
    mar = c(9, 3, 1, 1),
    mgp=c(1.5,0.4,0),
    tck=-0.008)
lapply(date.list,
       function(list.vector) {
         redox.plot(geochem.data.to.use = geochem.data,
                    date.range = list.vector,
                    RMs.to.use = c(286, 300, 310, 318),
                    nitrate.plotting.factor = nitrate.plotting.factor.to.use,
                    MeHg.plotting.factor = 2.5,
                    plot.Mn.instead.of.sulfide.YES.or.NO = "NO",
                    sulfide.plotting.factor = sulfide.plotting.factor.to.use,
                    max.depth = 80)
         
       })


#### 2018 profiles:  ####
all.dates[year(all.dates) == 2018]
date.list <- list(
  # c("2018-04-12", "2018-04-12"),
  # c("2018-05-01", "2018-05-01"),
  # c("2018-05-22", "2018-05-22"),
  c("2018-06-18", "2018-06-19"),
  # c("2018-07-10", "2018-07-10"),
  c("2018-07-31", "2018-07-31"),
  # c("2018-08-22", "2018-08-22"),
  # c("2018-09-11", "2018-09-11"),
  c("2018-09-24", "2018-09-26"),
  c("2018-10-16", "2018-10-16")
  # c("2018-11-06", "2018-11-06"),
  # c("2018-11-27", "2018-11-27"),
  # c("2018-12-11", "2018-12-11")
)
par(mfrow = c(4,4),
    mar = c(9, 3, 1, 1),
    mgp=c(1.5,0.4,0),
    tck=-0.008)
lapply(date.list,
       function(list.vector) {
         redox.plot(geochem.data.to.use = geochem.data,
                    date.range = list.vector,
                    RMs.to.use = c(286, 300, 310, 318),
                    nitrate.plotting.factor = nitrate.plotting.factor.to.use,
                    MeHg.plotting.factor = 2.5,
                    plot.Mn.instead.of.sulfide.YES.or.NO = "NO",
                    sulfide.plotting.factor = sulfide.plotting.factor.to.use,
                    max.depth = 80)
         
       })



dev.off()
