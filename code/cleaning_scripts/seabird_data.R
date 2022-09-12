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
  rename(depth = measurement_depth_m,
         elevation_m = brownlee_reservoir_measurement_elevation_m) %>%
  mutate(RM = as.character(snake_river_mile),
         date = as.Date(measurement_collection_date_mm_dd_yy_h_mm)) %>%
  filter(RM %in% RM.to.include) %>%
  select(RM, date, depth, elevation_m, water_temp_deg_c, diss_oxy_mg_per_l,
         spec_cond_ms_per_cm, ph, turb_ntu)
  

#### Read out data ####
saveRDS(seacat.data,
        "dataEdited/seabird/seabird_data.rds")
