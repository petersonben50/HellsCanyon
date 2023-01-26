#### code/geochem/HgT_time_course.R ####
# Benjamin D. Peterson


#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(ggpubr)
library(lubridate)
library(readxl)
library(tidyverse)
source("code/geochem/profile_functions.R")


#### Read in geochem data ####
geochem.data <- read.csv("dataEdited/geochem/geochem_WC.csv") %>%
  filter(depth != "SW") %>%
  mutate(depth = as.numeric(depth)) %>%
  arrange(depth)
geochem.data.of.interest <- geochem.data %>%
  mutate(RM_date = paste(RM, date, sep = "-")) %>%
  filter(RM_date %in% (geochem.data %>%
                         filter(month(date) >=5 & month(date) <= 11,
                                constituent == "HgT_diss_ngL") %>%
                         group_by(date, RM, constituent) %>%
                         summarise(count = n()) %>%
                         ungroup() %>%
                         filter(count >= 3) %>%
                         mutate(RM_date = paste(RM, date, sep = "-")) %>%
                         select(RM_date) %>%
                         unlist(use.names = FALSE)
                       )) %>%
  select(-RM_date) %>%
  filter(year(date) != 2015)
sampling.dates <- unique(geochem.data.of.interest$date)



#### Read in inflow data ####
inflow.data <- read_xlsx("dataRaw/dataRelease_chem/v2/HCC - V2 Data Release_Master File_Oct 2019 to present_V2_9_forModelingGroup.xlsx",
                              sheet = "Table_3_Water_V2") %>%
  # Ensure the date and time column is in correct format
  mutate(date = date(ymd_hms(sample_collection_date_mm_dd_yy_h_mm))) %>%
  # Rename some data columns
  rename(RM = snake_river_mile,
         depth = sample_depth_m,
         elevation_m = brownlee_reservoir_sample_elevation_m,
         HgT_diss_ngL = u_thg_ng_per_l) %>%
  # Select the data of interest
  select(RM, date, depth, elevation_m, HgT_diss_ngL) %>%
  gather(key = constituent,
         value = concentration,
         -c(1:4)) %>%
  filter(date %in% as.Date(sampling.dates),
         RM == 345.6) %>%
  mutate(elevation_m = max(geochem.data.of.interest$elevation_m) + 10,
         concentration = as.numeric(concentration),
         date = as.Date(date)) %>%
  filter(!is.na(concentration)) %>%
  select(date, elevation_m, constituent, concentration)



#### Generate needed color ramp ####
colorize.function = colorRampPalette(cb.translator[c("black", "orange", "yellow")])
nitrate.vector = colorize.function(7)
names(nitrate.vector) <- floor_date(seq(from = as.Date("2017-05-01"),
                                        to = as.Date("2017-11-01"),
                                        by = 1),
                                    'month') %>%
  unique() %>% month(label = TRUE)


#### Function to generate nitrate plots ####
function.of.function <- function(year.you.want) {
  time.course.profile.plot(geochem.data.to.use = geochem.data.of.interest,
                           years.to.use = year.you.want,
                           concentrations.to.use = c(-0.1, 7),
                           RMs.to.use = c(286,300, 310),
                           parameter.to.plot = "HgT_diss_ngL",
                           xlabel.to.use = "HgT (ng/L)",
                           color.vector.to.use = nitrate.vector,
                           inflow.data.to.use = inflow.data) +
    facet_wrap(~RM) +
    theme(text = element_text(size = 9))
}
plot.2016 <- function.of.function(2016) + theme(legend.position = "none")
plot.2017 <- function.of.function(2017) + theme(legend.position = "none")
plot.2018 <- function.of.function(2018) + theme(legend.position = "none")
plot.2019 <- function.of.function(2019) + theme(legend.position = "none")


#### Extract legend to plot separately ####
legend.to.use <- get_legend(function.of.function(2016))
# Convert to a ggplot and print
plotted.legend <- as_ggplot(legend.to.use)


#### Save out nitrate plots ####
pdf("results/geochem/HgT_time_course.pdf",
    width = 7,
    height = 8)
ggarrange(plot.2016, plot.2017, plot.2018, plot.2019,
          plotted.legend,
          labels = c("a", "b", "c", "d"),
          ncol = 2,
          nrow = 3)
dev.off()

