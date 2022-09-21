#### code/geochem/prepare_redox_grouping_data.R ####
# Benjamin D. Peterson

#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(tidyverse)
source("code/geochem/profile_functions.R")


#### Read in geochem data ####
geochem.data <- read.csv("dataEdited/geochem/geochem_WC.csv") %>%
  filter(constituent %in% c("MeHg_diss_ngL", "f_inorganic_sulfide_mg_per_l", "f_mn_mg_per_l",
                            "doc_boulder_mgc_per_l", "suva_254nm_l_per_mgc_per_m", "f_no3_mg_n_per_l")) %>%
  spread(key = constituent,
         value = concentration) %>%
  mutate(RM = as.character(RM),
         date = as.Date(date)) %>%
  # Data is super sparse before May, so let's not include
  # any of that
  filter(month(date) >= 5 &
           month(date) < 12)
# 968 observations.
# Keep only ones with dissolved MeHg measurements (hint: all of these do).
# Also only keep those with nitrate and Mn.
geochem.data <- geochem.data %>%
  filter(!is.na(MeHg_diss_ngL) &
           !is.na(f_no3_mg_n_per_l) &
           !is.na(f_mn_mg_per_l))
# 552 samples with nitrate measurements, of which 7 are
# missing Mn. These 7 include 2 that had measured (but
# non-detectable) sulfide, and none were high in MeHg.
geochem.data %>%
  filter(!is.na(f_inorganic_sulfide_mg_per_l)) %>%
  dim()
# Only 175 of the remaining samples had sulfide measured
geochem.data %>%
  filter(f_inorganic_sulfide_mg_per_l > 0.01) %>%
  dim()
# Of those 175 samples, only 4 had detectable sulfide (0.01 ppm).

# Don't keep the random ones we did from 305 and 314
geochem.data <- geochem.data %>%
  filter(RM != 314,
         RM != 305)


#### See where DO and geochem data matches up ####

# First, let's see where the data matches up.
seabird.data <- readRDS("dataEdited/seabird/seabird_data.rds") %>%
  filter(!is.na(diss_oxy_mg_per_l))

paired.list <- left_join(unique(geochem.data[, c("RM", "date")]),
                         unique(seabird.data[, c("RM", "date")]) %>%
                           mutate(seabird = "yes")) %>%
  arrange(seabird, date, RM)


#### Now find suitable replacements for DO measurements ####
seabird.data[(seabird.data$date == "2015-06-18" &
               seabird.data$RM == 300), "date"] <- as.Date("2015-06-04")
seabird.data[(seabird.data$date == "2015-06-18" &
               seabird.data$RM == 286), "date"] <- as.Date("2015-06-05")
seabird.data[(seabird.data$date == "2015-06-18" &
               seabird.data$RM == 310), "date"] <- as.Date("2015-06-06")
seabird.data[(seabird.data$date == "2015-09-23"), "date"] <- as.Date("2015-09-08")

# Used the profile from 2015-10-19 for both 2015-10-19 and the earlier trip
# (2015-09-29 to 2015-10-01)
duplicated.data <- seabird.data[(seabird.data$date == "2015-10-19"), ]
duplicated.data[(duplicated.data$date == "2015-10-19" &
                   duplicated.data$RM == 286), "date"] <- as.Date("2015-09-29")
duplicated.data[(duplicated.data$date == "2015-10-19" &
                   duplicated.data$RM == 300), "date"] <- as.Date("2015-09-30")
duplicated.data[(duplicated.data$date == "2015-10-19" &
                   duplicated.data$RM == 310), "date"] <- as.Date("2015-10-01")
seabird.data <- rbind(seabird.data,
                      duplicated.data)
rm(duplicated.data)
# No good analogue for 2015-11-16
geochem.data <- geochem.data %>% filter(date != "2015-11-16")

# 2016
seabird.data[(seabird.data$date == "2016-09-21" &
              seabird.data$RM == 286), "date"] <- as.Date("2016-10-03")
seabird.data[(seabird.data$date == "2016-09-21" &
               seabird.data$RM == 300), "date"] <- as.Date("2016-10-05")
seabird.data[(seabird.data$date == "2016-09-21" &
               seabird.data$RM %in% c(310, 318)), "date"] <- as.Date("2016-10-06")
seabird.data[(seabird.data$date == "2016-11-03" &
                seabird.data$RM == 286), "date"] <- as.Date("2016-11-15")
seabird.data[(seabird.data$date == "2016-11-03" &
                seabird.data$RM %in% c(300, 318)), "date"] <- as.Date("2016-11-22")

# 2017
geochem.data <- geochem.data %>% filter(!(RM == 318 & date == "2017-05-01"))
seabird.data[(seabird.data$date == "2017-06-14" &
               seabird.data$RM == "300"), "date"] <- as.Date("2017-06-06")
geochem.data <- geochem.data %>% filter(date != "2017-11-14")

# 2018
seabird.data[(seabird.data$date == "2018-06-12" &
               seabird.data$RM %in% c(286, 300)), "date"] <- as.Date("2018-06-18")
seabird.data[(seabird.data$date == "2018-06-12" &
               seabird.data$RM %in% c(310, 318)), "date"] <- as.Date("2018-06-19")
geochem.data <- geochem.data %>% filter(date != "2018-09-26")
geochem.data

seabird.data[(seabird.data$date == "2019-07-23" &
               seabird.data$RM == 286), "date"] <- as.Date("2019-07-24")
seabird.data[(seabird.data$date == "2019-07-18" &
               seabird.data$RM == 300), "date"] <- as.Date("2019-07-25")


#### Scripts useful for testing above conversions ####
# par(mfrow = c(1,4),
#     mar = c(9, 4.5, 1, 1),
#     mgp=c(1.5,0.4,0),
#     tck=-0.008)
# redox.plot(geochem.data.to.use = geochem.data %>% gather(key = constituent, value = concentration, -c(1:4)),
#            date.range = c("2018-06-17", "2018-06-23"),
#            DO.data.to.use = seabird.data,
#            DO.date = "2018-06-12",
#            RMs.to.use = c(286, 300, 310, 318),
#            nitrate.plotting.factor = 5,
#            MeHg.plotting.factor = 2.5,
#            sulfide.plotting.factor = 10,
#            plot.Mn.instead.of.sulfide.YES.or.NO = "NO",
#            plot.Mn.with.sulfide = "YES",
#            Mn.plotting.factor = 5,
#            plot.DO = "YES")
# left_join(unique(geochem.data[, c("RM", "date")]),
#                          unique(seabird.data[, c("RM", "date")]) %>%
#                            mutate(seabird = "yes")) %>%
#   arrange(seabird, date, RM) %>%
#   filter(is.na(seabird))
# par(mfrow = c(1,1))
# plot(log(geochem.data$f_mn_mg_per_l, 10), log(geochem.data$MeHg_diss_ngL, 10))


# Next we need to interpolate the data to match the elevations
# where we collected samples.
unique.dates.sites <- unique(geochem.data[, c("RM", "date")])
for (entry in 1:dim(unique.dates.sites)[1]) {
  seabird.data.temporary <- seabird.data %>%
    filter(RM == unique.dates.sites[entry, "RM"] &
             date == unique.dates.sites[entry, "date"])
  
  elevations.needed <- geochem.data %>%
    filter(RM == unique.dates.sites[entry, "RM"] &
             date == unique.dates.sites[entry, "date"]) %>%
    ungroup() %>%
    select(depth, elevation_m) %>%
    unique()
  
  interpol.data.DO <- approx(x = seabird.data.temporary$elevation_m,
                             y = seabird.data.temporary$diss_oxy_mg_per_l,
                             xout = as.numeric(elevations.needed$elevation_m),
                             rule = 2)
  seabird.data.adjusted.temp <- data.frame(RM = unique.dates.sites[entry, "RM"],
                                          date = unique.dates.sites[entry, "date"],
                                          depth = as.numeric(elevations.needed$depth),
                                          elevation_m = interpol.data.DO$x,
                                          diss_oxy_mg_per_l = interpol.data.DO$y)
  if (entry == 1) {
    seabird.data.adjusted <- seabird.data.adjusted.temp
  } else {
    seabird.data.adjusted <- rbind(seabird.data.adjusted,
                                  seabird.data.adjusted.temp)
  }
  rm(seabird.data.adjusted.temp, interpol.data.DO, elevations.needed, seabird.data.temporary)
}



#### Combine all data ####
all.redox.data <- full_join(geochem.data,
                            seabird.data.adjusted)


#### Assign redox state ####
all.redox.data <- all.redox.data %>%
  mutate(redox_status = "none_assigned")

# All samples where sulfide is detected are labeled as sulfidic
all.redox.data[which(all.redox.data$f_inorganic_sulfide_mg_per_l > 0.01), "redox_status"] <- "sulfidic"
# 
all.redox.data[which(all.redox.data$diss_oxy_mg_per_l > 1 &
                       all.redox.data$f_mn_mg_per_l < 0.1), "redox_status"] <- "oxic"
# If nitrate is below 0.05 (the max DDL), there's no DO, and
# sulfide was non-detect, we'll call that no nitrate no sulfide
all.redox.data[which(all.redox.data$f_no3_mg_n_per_l < 0.05 &
                       all.redox.data$diss_oxy_mg_per_l < 1 &
                       all.redox.data$f_inorganic_sulfide_mg_per_l == 0.01), "redox_status"] <- "no_nitrate_no_sulfide"
# If nitrate is below 0.05 (the max DDL), there's no DO, and
# sulfide was non-detect, we'll call that no nitrate possible sulfide
all.redox.data[which(all.redox.data$f_no3_mg_n_per_l < 0.05 &
                       all.redox.data$diss_oxy_mg_per_l < 1 &
                       is.na(all.redox.data$f_inorganic_sulfide_mg_per_l)), "redox_status"] <- "no_nitrate_possible_sulfide"
# All the rest are suboxic
all.redox.data[which(all.redox.data$redox_status == "none_assigned"), "redox_status"] <- "suboxic"


#### Save out data ####
saveRDS(all.redox.data,
        "dataEdited/geochem/geochem_WC_adjusted_for_redox_classification.rds")
