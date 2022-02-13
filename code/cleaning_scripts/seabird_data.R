#### code/cleaning_scripts/seabird_data.R ####
# Benjamin D. Peterson


#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(readxl)
library(tidyverse)


#### Read in data ####
# 2017 and 2018 data
sonde.data.2017.2018 <- read_xlsx("dataRaw/SeaCat/sondeProfiles_intensives2017-2018.xlsx") %>%
  select(-c(`ChlA (EXO)`, `Turbidity (EXO)`, `FOM(EXO)`))
names(sonde.data.2017.2018) <- c("date", "RM", "depth", "spec_cond", "redox_pot", "pH",
                                 "temp", "DO", "ChlA", "turb", "PAR", "cDOM")
sonde.data.2017.2018 <- sonde.data.2017.2018 %>%
  mutate(date = as.Date(date)) %>%
  filter(month(date) == 9)
# 2019 data
sonde.data.2019 <- read_xlsx("dataRaw/SeaCat/2019_intensive_seacat_brownlee.xlsx") %>%
  select(-c(Datetime, Time, Location))
names(sonde.data.2019) <- c("date", "depth", "temp", "pH", "redox_pot", "cDOM",
                            "ChlA", "turb", "PAR", "spec_cond", "DO", "RM")
sonde.data.2019 <- sonde.data.2019 %>%
  filter(RM != "322") %>%
  mutate(date = mdy(date)) %>%
  select(c("date", "RM", "depth", "spec_cond", "redox_pot", "pH",
           "temp", "DO", "ChlA", "turb", "PAR", "cDOM"))


#### Combine data ####
sonde.data <- rbind(sonde.data.2017.2018,
                    sonde.data.2019) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:3)) %>%
  select(RM, date, depth, constituent, concentration)
rm(sonde.data.2017.2018,
   sonde.data.2019)


#### Read out data ####
saveRDS(sonde.data,
        "dataEdited/seabird/seabird_data.rds")
