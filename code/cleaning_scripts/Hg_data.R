#### code/cleaning_scripts/water_chem_data_cleaning_2015_2018.R ####
# Benjamin D. Peterson

# This file contains the code to take the Hg data
# that Brett sent me from 2017, 2018, and 2019 and
# clean it up.


#### Get a haircut, and get a real job ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(dplyr)
library(lubridate)
library(tidyr)
library(readxl)



#### Read in 2015-2018 data ####
chem.286.data <- read_xlsx("dataRaw/waterChemistry/Poulin_Hells Canyon Data for Ben Peterson_07082019.xlsx",
                           sheet = 2) %>%
  sapply(as.character) %>%
  as.data.frame(stringsAsFactors = FALSE)
chem.300.data <- read_xlsx("dataRaw/waterChemistry/Poulin_Hells Canyon Data for Ben Peterson_07082019.xlsx",
                           sheet = 3) %>%
  sapply(as.character) %>%
  as.data.frame(stringsAsFactors = FALSE)
chem.310.data <- read_xlsx("dataRaw/waterChemistry/Poulin_Hells Canyon Data for Ben Peterson_07082019.xlsx",
                           sheet = 4) %>%
  sapply(as.character) %>%
  as.data.frame(stringsAsFactors = FALSE)
chem.318.data <- read_xlsx("dataRaw/waterChemistry/Poulin_Hells Canyon Data for Ben Peterson_07082019.xlsx",
                           sheet = 5) %>%
  sapply(as.character) %>%
  as.data.frame(stringsAsFactors = FALSE)


#### Clean 2015-2018 data ####
clean.chem.data <- chem.286.data %>%
  bind_rows(chem.300.data) %>%
  bind_rows(chem.310.data) %>%
  bind_rows(chem.318.data) %>%
  filter(!is.na(Site))
names(clean.chem.data) <- c("RM", "date", "depth", "TSS_mg.L", "DOC_mC.L",
                            "SUVA254_L.mg.m", "FTHG", "FMHG", "per_FMeHg",
                            "PTHG", "PMHG", "per_PMeHg", "Fe_ug.L",
                            "Fe.part_ug.L", "Mn_ug.L", "Mn.part_ug.L", "SO4_mg.L",
                            "NO3_mgN.L", "NH4_mg.L", "Cl_mg.L", "SO4_div_Cl_mg.mg",
                            "NO3_div_Cl_mg.mg", "sulfide_mg.L", "thiosulfate_mg.L")

clean.chem.data <- clean.chem.data %>%
  mutate(date = round(as.Date(date)))


#### Keep only Hg data from 2017 and 2018 intensives ####
Hg.2017.2018 <- clean.chem.data %>%
  select(RM, depth, date, FTHG, FMHG,
         PTHG, PMHG) %>%
  filter(year(date) %in% c(2017, 2018),
         month(date) == 9,
         day(date) > 15)
rm(list = ls(pattern = "chem"))




#### Read in and process 2019 data ####
Hg.data <- read_xlsx("dataRaw/waterChemistry/Hg/HCC_July 2019 water column.xlsx",
                     sheet = "results_2020-01-29")
Hg.2019 <- Hg.data %>%
  filter(depth > 0) %>%
  filter(medium %in% c("WS", "WSQ")) %>%
  mutate(RM = gsub("SR RM ", "", site_name))%>%
  select(RM, depth, sample_date, medium, constituent, result_value) %>%
  mutate(date = as.Date(sample_date)) %>%
  spread(key = constituent,
         value = result_value) %>%
  select(RM, depth, date, FTHG, 
         FMHG, PTHG, PMHG)
rm(Hg.data)


#### Combine data ####
all.Hg.data <- rbind(Hg.2017.2018,
                     Hg.2019)
rm(Hg.2017.2018,
   Hg.2019)


#### Calculate inorganic Hg and percent MeHg ####
all.Hg.data <- all.Hg.data %>%
  mutate(PiHg = as.numeric(PTHG) - as.numeric(PMHG),
         FiHg = as.numeric(FTHG) - as.numeric(FMHG)) %>%
  mutate(PMHG_per = as.numeric(PMHG) / as.numeric(PTHG),
         FMHG_per = as.numeric(FMHG) / as.numeric(FTHG))


#### Save out data ####
write.csv(all.Hg.data,
          file = "dataEdited/waterChemistry/Hg_data.csv",
          row.names = FALSE)
