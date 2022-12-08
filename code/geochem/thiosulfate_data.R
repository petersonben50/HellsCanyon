#### code/geochem/thiosulfate_data.R ####
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
thiosulfate.data <- read_xlsx("dataRaw/dataRelease_chem/v2/HCC - V2 Data Release_Master File_Oct 2019 to present_V2_9_forModelingGroup.xlsx",
                              sheet = "Table_3_Water_V2") %>%
  # Ensure the date and time column is in correct format
  mutate(date = date(ymd_hms(sample_collection_date_mm_dd_yy_h_mm))) %>%
  rename(RM = snake_river_mile,
         depth = sample_depth_m,
         elevation_m = brownlee_reservoir_sample_elevation_m,
         f_s2o3_mg_per_l = f_s2o3_mg_per_l_) %>%
  select(date, RM, depth, elevation_m,
         f_s2o3_mg_per_l, f_s2o3_mg_per_l_qa) %>%
  filter(date > "2016-01-01",
         RM %in% c(286, 300, 310),
         depth != "--",
         f_s2o3_mg_per_l != "--") %>%
  mutate(depth = as.numeric(depth),
         elevation_m = as.numeric(elevation_m))



#### Read in sulfide and nitrate data ####
geochem.data <- read.csv("dataEdited/geochem/geochem_WC.csv") %>%
  filter(constituent %in% c("f_inorganic_sulfide_mg_per_l",
                            "f_no3_mg_n_per_l")) %>%
  spread(key = constituent,
         value = concentration) %>%
  mutate(date = as.Date(date),
         RM = as.character(RM))


#### Combine data ####
S.data <- left_join(thiosulfate.data, geochem.data)

S.data.thiosulfate.detected <- S.data %>%
  filter(!grepl("<", f_s2o3_mg_per_l,)) %>%
  mutate(f_s2o3_mg_per_l = round(as.numeric(f_s2o3_mg_per_l), 2)) %>%
  arrange(date, RM, depth)
S.data.thiosulfate.detected

