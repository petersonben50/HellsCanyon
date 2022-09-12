#### code/cleaning_scripts/water_chem_data.R ####
# Benjamin D. Peterson

# This file contains the code to pull out the data I
# need from the data release and convert it into the
# form I need.
# It also takes the nitrite data that Brett sent me
# and adds that to the other data.


#### Get a haircut, and get a real job ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(dplyr)
library(lubridate)
library(tidyr)
library(readxl)


#### Read in data release data ####
all.geochem.data <- read_xlsx("dataRaw/dataRelease_chem/v2/HCC - V2 Data Release_Master File_Oct 2019 to present_V2_9_forModelingGroup.xlsx",
                              sheet = "Table_3_Water_V2") %>%
  # Ensure the date and time column is in correct format
  mutate(date = date(ymd_hms(sample_collection_date_mm_dd_yy_h_mm))) %>%
  # Rename some data columns
  rename(RM = snake_river_mile,
         depth = sample_depth_m,
         elevation_m = brownlee_reservoir_sample_elevation_m,
         HgT_diss_ngL = f_thg_ng_per_l,
         MeHg_diss_ngL = f_mehg_ng_per_l,
         MeHg_diss_percent = percent_f_mehg,
         HgT_part_ngL = p_thg_vol_ng_per_l,
         MeHg_part_ngL = p_mehg_vol_ng_per_l,
         MeHg_part_percent = percent_p_mehg) %>%
  mutate(f_mn_mg_per_l = (as.numeric(f_mn_mcg_per_l) / 1000)) %>%
  # Select the data of interest
  select(RM, date, depth, elevation_m, medium_code,
         HgT_diss_ngL, MeHg_diss_ngL, MeHg_diss_percent,
         HgT_part_ngL, MeHg_part_ngL, MeHg_part_percent,
         f_mn_mg_per_l, f_so4_mg_per_l, f_no3_mg_n_per_l,
         f_cl_mg_per_l, f_inorganic_sulfide_mg_per_l) %>%
  # Make it long!
  gather(key = constituent,
         value = concentration,
         -c(1:5)) %>%
  filter(depth != "--",
         depth != "SW",
         elevation_m != "--",
         RM != "--",
         concentration != "--") %>%
  mutate(concentration = gsub("<", "", concentration)) %>%
  mutate(concentration = as.numeric(concentration),
         depth = as.numeric(depth))


#### Add nitrite data for 2019 ####
nitrite.2019.data <- read_xlsx("dataRaw/waterChemistry/HCC_01272020_July 2019 Intensive Metals Data_For Ben Peterson.xlsx",
                                  skip = 2) %>%
  rename(depth = Depth) %>%
  filter(!is.na(RM),
         !is.na(NO2)) %>%
  mutate(date = as.Date(as.numeric(SampDate),
                        origin = "1899-12-30")) %>%
  mutate(medium_code = "WS") %>%
  select(RM, date, depth, medium_code, NO2) %>%
  mutate(NO2 = gsub(pattern = "na", replacement = 0,
                    NO2),
         NO2 = gsub(pattern = "-2E-3", replacement = 0,
                    NO2),
         NO2 = as.numeric(NO2)) %>%
  gather(key = constituent,
         value = concentration,
         -c(RM, depth, date, medium_code)) %>%
  filter(depth != "na",
         !grepl("cm", depth)) %>%
  mutate(depth = as.numeric(depth)) %>%
  left_join(all.geochem.data %>%
              select(RM, depth, date, elevation_m)) %>%
  mutate(elevation_m = as.character(as.numeric(elevation_m) %>% round(digits = 1))) %>%
  select(RM, date, depth, elevation_m, medium_code,
         constituent, concentration)



#### Calculate inorganic concentration ####
inorganic.Hg.data <- all.geochem.data %>%
  filter(constituent %in% c("HgT_diss_ngL", "HgT_part_ngL", "MeHg_diss_ngL", "MeHg_part_ngL")) %>%
  group_by(RM, depth, date, elevation_m, medium_code, constituent) %>%
  summarise(concentration = mean(concentration)) %>%
  ungroup() %>%
  spread(key = constituent,
         value = concentration) %>%
  mutate(iHg_diss_ngL = HgT_diss_ngL - MeHg_diss_ngL,
         iHg_part_ngL = HgT_part_ngL - MeHg_part_ngL) %>%
  gather(key = constituent,
         value = concentration,
         -c(RM, depth, date, elevation_m, medium_code)) %>%
  filter(constituent %in% c("iHg_diss_ngL", "iHg_part_ngL")) %>%
  filter(!is.na(concentration))



#### Combine the nitrite data with the other data
all.geochem.data <- rbind(all.geochem.data,
                          nitrite.2019.data,
                          inorganic.Hg.data) %>%
  group_by(RM, date, depth, elevation_m, constituent) %>%
  summarise(concentration = mean(concentration)) %>%
  ungroup() %>%
  mutate(elevation_m = as.numeric(elevation_m))





#### Add DO data ####
# First need to interpolate the data to match the elevations
# where we collected samples.
seacat.data <- readRDS("dataEdited/seabird/seabird_data.rds") %>%
  filter(!is.na(diss_oxy_mg_per_l))

unique.dates.sites.geochem <- unique(all.geochem.data[, c("RM", "date")]) %>%
  as.data.frame()
unique.dates.sites.seacat <- unique(seacat.data[, c("RM", "date")]) %>%
  as.data.frame()
unique.dates.sites <- inner_join(unique.dates.sites.geochem,
                                 unique.dates.sites.seacat) %>%
  arrange(date)

for (entry in 1:dim(unique.dates.sites)[1]) {
  seacat.data.temporary <- seacat.data %>%
    filter(RM == unique.dates.sites[entry, "RM"] &
             date == unique.dates.sites[entry, "date"])
  
  elevations.needed <- all.geochem.data %>%
    filter(RM == unique.dates.sites[entry, "RM"] &
             date == unique.dates.sites[entry, "date"]) %>%
    ungroup() %>%
    select(depth, elevation_m) %>%
    unique()
  
  interpol.data.DO <- approx(x = seacat.data.temporary$elevation_m,
                             y = seacat.data.temporary$diss_oxy_mg_per_l,
                             xout = as.numeric(elevations.needed$elevation_m),
                             rule = 2)
  seacat.data.adjusted.temp <- data.frame(RM = unique.dates.sites[entry, "RM"],
                                          date = unique.dates.sites[entry, "date"],
                                          depth = as.numeric(elevations.needed$depth),
                                          elevation_m = interpol.data.DO$x,
                                          diss_oxy_mg_per_l = interpol.data.DO$y)
  if (entry == 1) {
    seacat.data.adjusted <- seacat.data.adjusted.temp
  } else {
    seacat.data.adjusted <- rbind(seacat.data.adjusted,
                                  seacat.data.adjusted.temp)
  }
}

seacat.data.adjusted <- seacat.data.adjusted %>%
  rename(concentration = diss_oxy_mg_per_l) %>%
  mutate(constituent = "diss_oxy_mg_per_l") %>%
  select(RM, date, depth, elevation_m, constituent, concentration)

all.geochem.data.final <- rbind(all.geochem.data,
                                seacat.data.adjusted)


#### Write out water column dissolved data ####
write.csv(all.geochem.data.final,
          "dataEdited/waterChemistry/geochem_WC.csv",
          row.names = FALSE,
          quote = FALSE)
