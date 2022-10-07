#### code/geochem/calculate_age_anoxia.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Read in geochem data ####
# Only use the data that was used in Mn vs. MeHg figure.
all.redox.data <- readRDS("dataEdited/geochem/geochem_WC_adjusted_for_redox_classification.rds") %>%
  mutate(sampling.year = as.character(year(date)))
# Identify unique depths at each RM within a year
unique.locations.per.year <- all.redox.data %>%
  group_by(sampling.year, RM, depth) %>%
  summarize(times.sampled = n()) %>%
  ungroup() %>%
  select(sampling.year, RM, depth, times.sampled)


#### See where DO and geochem data matches up ####

# First, let's see where the data matches up.
seabird.data <- readRDS("dataEdited/seabird/seabird_data.rds") %>%
  filter(!is.na(diss_oxy_mg_per_l)) %>%
  select(RM, date, depth, diss_oxy_mg_per_l)



#### Estimate initial time of anoxia at each of the sampling locations ####
for (entry in 1:dim(unique.locations.per.year)[1]) {

  sampling.year.of.interest <- as.character(unique.locations.per.year[entry, "sampling.year"])
  RM.of.interest <- as.character(unique.locations.per.year[entry, "RM"])
  depth.of.interest <- unique.locations.per.year[entry, "depth"] %>% unlist(use.names = FALSE)

  data.subset <- seabird.data %>%
    filter(year(date) == sampling.year.of.interest,
           RM == RM.of.interest)

  for (date in unique(data.subset$date)) {

    list.of.dates <- unique(data.subset$date)
    DO.at.depth <- sapply(list.of.dates,
                          function(date.to.use) {
                            subset.of.data.subset <- data.subset %>%
                              filter(date == date.to.use)
                            approx(x = subset.of.data.subset$depth,
                                   y = subset.of.data.subset$diss_oxy_mg_per_l,
                                   xout = depth.of.interest,
                                   rule = 2)$y
                          })
    DO.at.depth.by.date <- data.frame(date = list.of.dates,
                                      DO = DO.at.depth)

  }

  interp.DO.by.date <- approx(x = DO.at.depth.by.date$date,
                              y = DO.at.depth.by.date$DO,
                              xout = seq(min(DO.at.depth.by.date$date),
                                         max(DO.at.depth.by.date$date),
                                         by = 1),
                              rule = 2)
  date.of.initial.anoxia.df <- data.frame(date = as.Date(interp.DO.by.date$x),
                                       DO = interp.DO.by.date$y) %>%
    mutate(anoxic = (DO < 1)) %>%
    group_by(anoxic) %>%
    summarise(date = min(date)) %>%
    ungroup()
  if (any(date.of.initial.anoxia.df$anoxic)) {
    date.of.initial.anoxia.df <- date.of.initial.anoxia.df %>%
      filter(anoxic == TRUE) %>%
      select(date) %>%
      mutate(RM = RM.of.interest,
             sampling.year = sampling.year.of.interest,
             depth = depth.of.interest,
             date.of.initial.anoxia = date) %>%
      select(RM, depth, sampling.year, date.of.initial.anoxia)
  } else {
    date.of.initial.anoxia.df <- data.frame(RM = RM.of.interest,
                                            depth = depth.of.interest,
                                            sampling.year = sampling.year.of.interest,
                                            date.of.initial.anoxia = as.Date(paste(sampling.year.of.interest, "-12-31",
                                                                                   sep = ""))
                                            )
  }

  if (entry == 1) {
    anoxia.start.df <- date.of.initial.anoxia.df
  } else {
    anoxia.start.df <- rbind(anoxia.start.df,
                             date.of.initial.anoxia.df)
  }
}
saveRDS(anoxia.start.df,
        "dataEdited/geochem/anoxia_start_dates.rds")

