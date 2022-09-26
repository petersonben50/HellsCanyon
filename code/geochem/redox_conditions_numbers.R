#### code/geochem/redox_conditions_numbers.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
# library(lubridate)
library(tidyverse)
# cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Read in data ####
geochem.data.adj <- readRDS("dataEdited/geochem/geochem_WC_adjusted_for_redox_classification.rds")


# Total samples
dim(geochem.data.adj)[1]


# Total oxygen-depleted samples
geochem.data.adj %>%
  filter(redox_status != "oxic") %>%
  dim()


# Total suboxic samples
geochem.data.adj %>%
  filter(redox_status == "suboxic") %>%
  dim()


# Total number of samples with no nitrate
geochem.data.adj %>%
  filter(redox_status %in% c("no_nitrate_no_sulfide",
                             "no_nitrate_possible_sulfide",
                             "sulfidic")) %>%
  dim()
  

# Total number of samples with no sulfide
geochem.data.adj %>%
  filter(redox_status %in% c("sulfidic")) %>%
  dim()
