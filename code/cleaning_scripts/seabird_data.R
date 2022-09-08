#### code/cleaning_scripts/seabird_data.R ####
# Benjamin D. Peterson


#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(readxl)
library(tidyverse)


#### Define river miles of interest ####
RM.to.include <- c(286, 300, 310, 318)

#### Read in data ####
seacat.data <- read_xlsx("dataRaw/dataRelease_chem/v2/HCC - V2 Data Release_Master File_Oct 2019 to present_V2_9_forModelingGroup.xlsx",
                         sheet = "Table_6_Brownlee_Profiles") %>%
  rename(RM = snake_river_mile,
         elevation_m = brownlee_reservoir_measurement_elevation_m) %>%
  mutate(date = as.Date(measurement_collection_date_mm_dd_yy_h_mm)) %>%
  filter(RM %in% RM.to.include) %>%
  select(RM, date, elevation_m, diss_oxy_mg_per_l) %>%
  filter(!is.na(diss_oxy_mg_per_l))

unique.dates.sites <- unique(seacat.data[, c("RM", "date")]) %>%
  as.data.frame()
  

for (entry in 1:dim(unique.dates.sites)[1]) {
  seacat.data.temp <- seacat.data %>%
    filter(RM == unique.dates.sites[entry, "RM"] &
             date == unique.dates.sites[entry, "date"])
  interpol.data.DO <- approx(x = seacat.data.temp$elevation_m,
                             y = seacat.data.temp$diss_oxy_mg_per_l,
                             xout = seq(floor(min(seacat.data.temp$elevation_m)),
                                        ceiling(max(seacat.data.temp$elevation_m)),
                                        1),
                             rule = 2)
  seacat.data.adjusted.temp <- data.frame(RM = unique.dates.sites[entry, "RM"],
                                          date = unique.dates.sites[entry, "date"],
                                          elevation_m = interpol.data.DO$x,
                                          diss_oxy_mg_per_l = interpol.data.DO$y)
  if (entry == 1) {
    seacat.data.adjusted <- seacat.data.adjusted.temp
  } else {
    seacat.data.adjusted <- rbind(seacat.data.adjusted,
                                  seacat.data.adjusted.temp)
  }
}


#### Read out data ####
saveRDS(seacat.data.adjusted,
        "dataEdited/seabird/seabird_data.rds")
