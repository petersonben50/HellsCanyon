#### code/cleaning_scripts/water_chem_data.R ####
# Benjamin D. Peterson

# This file contains the code to take the geochem data
# that Brett sents me and convert it into an easily 
# parseable format. 


#### Get a haircut, and get a real job ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(dplyr)
library(lubridate)
library(tidyr)
library(readxl)


#### Read in data from 2018 and earlier ####
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


#### Clean data from 2018 and earlier ####
clean.chem.data <- chem.286.data %>%
  bind_rows(chem.300.data) %>%
  bind_rows(chem.310.data) %>%
  bind_rows(chem.318.data) %>%
  filter(!is.na(Site))
rm(chem.286.data, chem.300.data,
   chem.310.data, chem.318.data)
names(clean.chem.data) <- c("RM", "date", "depth", "TSS_mg.L", "DOC_mC.L",
                            "SUVA254_L.mg.m", "FTHg_ng.L", "FMeHg_ng.L", "per_FMeHg",
                            "PTHg_ng.L", "PMeHg_ng.L", "per_PMeHg", "Fe_ug.L",
                            "Fe.part_ug.L", "Mn_ug.L", "Mn.part_ug.L", "SO4_mg.L",
                            "NO3_mgN.L", "NH4_mg.L", "Cl_mg.L", "SO4_div_Cl_mg.mg",
                            "NO3_div_Cl_mg.mg", "sulfide_mg.L", "thiosulfate_mg.L")
# Replace the text strings in the Mn particulate with 0
clean.chem.data$Mn.part_ug.L[grep("<", clean.chem.data$Mn.part_ug.L)] <- NA

clean.chem.data <- clean.chem.data %>%
  mutate(date = round(as.Date(date)))


#### Replace nds with 0s ####
clean.chem.data$sulfide_mg.L[which(clean.chem.data$sulfide_mg.L == "nd")] <- 0


#### Clean the data ####
clean.chem.data.201X <- clean.chem.data %>%
  select(RM, depth, date, Fe.part_ug.L, Fe_ug.L,
         Mn.part_ug.L, Mn_ug.L, NO3_mgN.L, NH4_mg.L,
         SO4_mg.L, sulfide_mg.L, Cl_mg.L) %>%
  gather(key = constituent,
         value = concentration,
         -c(RM, depth, date)) %>%
  filter(!is.na(concentration))
rm(clean.chem.data)









####################################################################
#### Metals in water column data from 2019 intensive from Brett ####

#### Read in data ####
# Huge data file incoming! Let's check it out. 
# Could probably just use the nitrate/nitrite/sulfate data from this as well...
water.chem.data.2019 <- read_xlsx("dataRaw/waterChemistry/HCC_01272020_July 2019 Intensive Metals Data_For Ben Peterson.xlsx",
                                  skip = 2)

#### Clean up data ####
# Separate out dissolved and particulate.
# We also don't need to include the inflow grab
# samples (I think that's what they are).
# Also separate out the WC from the OL samples
water.chem.data.dissolved.WC <- water.chem.data.2019[3:44, ] %>%
  mutate(date = as.Date(as.numeric(SampDate), origin = "1899-12-30"),
         depth = Depth,
         Fe_ug.L = round(as.numeric(Fe), 3),
         Mn_ug.L = round(as.numeric(Mn), 3),
         S_mg.L = round(as.numeric(S), 3),
         Cl_mg.L = round(as.numeric(Cl), 3),
         SO4_mg.L = round(as.numeric(SO4), 3),
         NO3_mgN.L = round(as.numeric(NO3), 3),
         NO2_mgN.L = round(as.numeric(NO2), 3),
         PO4_mgP.L = round(as.numeric(PO4), 3))
water.chem.data.dissolved.WC$sulfide_mg.L <- rep(NA,
                                                 length(water.chem.data.dissolved.WC$FieldID))
water.chem.data.dissolved.WC <- water.chem.data.dissolved.WC %>%
  select(FieldID, date, RM, depth,
         Fe_ug.L, Mn_ug.L,S_mg.L,
         Cl_mg.L, SO4_mg.L,sulfide_mg.L,
         NO3_mgN.L, NO2_mgN.L, PO4_mgP.L)

#### Fill in the depth for SnakeR-RM 310-26-REP m ####
water.chem.data.dissolved.WC$depth[which(water.chem.data.dissolved.WC$FieldID == "SnakeR-RM 310-26-REP m")] <- 26
water.chem.data.dissolved.WC$RM[which(water.chem.data.dissolved.WC$FieldID == "SnakeR-RM 310-26-REP m")] <- 310




#### Pull out the whole water data ####
water.chem.data.2019.whole <- water.chem.data.2019[77:118, ] %>%
  filter(FieldID != "na") %>%
  mutate(date = as.Date(as.numeric(SampDate), origin = "1899-12-30"),
         depth = Depth,
         Fe.whole_ug.L = round(as.numeric(Fe), 3),
         Mn.whole_ug.L = round(as.numeric(Mn), 3),
         Cl.whole_mg.L = round(as.numeric(Cl), 3),
         S.whole_mg.L = round(as.numeric(S), 3)) %>%
  select(FieldID, date, RM, depth, Fe.whole_ug.L, Mn.whole_ug.L, S.whole_mg.L, Cl.whole_mg.L) 


#### Combine data and calculate the particulate fraction ####
water.chem.data.2019.all <- full_join(water.chem.data.2019.whole,
                                      water.chem.data.dissolved.WC) %>%
  filter(RM != "na") %>%
  mutate(Fe.part_ug.L = Fe.whole_ug.L - Fe_ug.L,
         Mn.part_ug.L = Mn.whole_ug.L - Mn_ug.L,
         Cl.part_mg.L = Cl.whole_mg.L - Cl_mg.L,
         S.part_mg.L = S.whole_mg.L - S_mg.L)
rm(water.chem.data.2019.whole,
   water.chem.data.dissolved.WC,
   water.chem.data.2019)


#### Final set of 2019 data ####
clean.chem.data.2019 <- water.chem.data.2019.all %>%
  select(RM, depth, date,
         Fe_ug.L, Fe.part_ug.L, Mn_ug.L, Mn.part_ug.L,
         S_mg.L, Cl_mg.L, SO4_mg.L, sulfide_mg.L,
         NO3_mgN.L, NO2_mgN.L, PO4_mgP.L) %>%
  gather(key = constituent,
         value = concentration,
         -c(RM, depth, date)) %>%
  filter(!is.na(concentration))
rm(water.chem.data.2019.all)


#### Combine and clean data ####
clean.chem.data <- rbind(clean.chem.data.201X,
                         clean.chem.data.2019)
rm(clean.chem.data.201X,
   clean.chem.data.2019)




#### Write out water column dissolved data ####
write.csv(clean.chem.data,
          "dataEdited/waterChemistry/geochem_WC.csv",
          row.names = FALSE,
          quote = FALSE)
