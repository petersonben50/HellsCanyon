#### code/geochem/corewater.R ####
# Written by Benjamin D. Peterson


#### Get ready. Get set ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(ggpubr)
library(lubridate)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")

#### Read in corewater data ####
corewater.data <- read_xlsx("dataRaw/dataRelease_chem/v2/HCC - V2 Data Release_Master File_Oct 2019 to present_V2_9_forModelingGroup.xlsx",
                            sheet = "Table_3_Water_V2") %>%
  filter(min_sample_height_above_sediment_m != "--",
         snake_river_mile != "--") %>%
  mutate(RM = as.character(snake_river_mile),
         date = as.Date(sample_collection_date_mm_dd_yy_h_mm),
         upper_height = as.numeric(max_sample_height_above_sediment_m),
         lower_height = as.numeric(min_sample_height_above_sediment_m),
         FTHG = as.numeric(f_thg_ng_per_l),
         FMHG = as.numeric(f_mehg_ng_per_l)) %>%
  mutate(height = (upper_height + lower_height) / 2) %>%
  select(RM, date, height,
         FTHG, FMHG) %>%
  filter(RM != 248.1)


#### Prepare water column data for sample immediately overlying corewater ####
unique.sample.locations <- corewater.data %>%
  select(RM, date) %>%
  unique() %>%
  mutate(corewater = "yes")
water.column.data <- read_xlsx("dataRaw/dataRelease_chem/v2/HCC - V2 Data Release_Master File_Oct 2019 to present_V2_9_forModelingGroup.xlsx",
                               sheet = "Table_3_Water_V2") %>%
  mutate(RM = as.character(snake_river_mile),
         date = as.Date(sample_collection_date_mm_dd_yy_h_mm),
         depth = sample_depth_m,
         elevation_m = brownlee_reservoir_sample_elevation_m,
         FTHG = as.numeric(f_thg_ng_per_l),
         FMHG = as.numeric(f_mehg_ng_per_l)) %>%
  filter(RM != "--",
         depth != "--") %>%
  select(RM, date, depth, elevation_m, FTHG, FMHG) %>%
  mutate(elevation_m = as.numeric(elevation_m)) #%>%
# Corewater at 305 from July 2019 intensive was on 2019-07-24, but profile was on 2019-07-22
water.column.data[which(water.column.data$RM == 305 &
                          water.column.data$date == "2019-07-22"), "date"] <- as.Date("2019-07-24")
# Corewater at 314 from July 2019 intensive was on 2019-07-23, but profile was on 2019-07-22
water.column.data[which(water.column.data$RM == 314 &
                          water.column.data$date == "2019-07-22"), "date"] <- as.Date("2019-07-23")
# Corewater at 318 from July 2019 intensive was on 2019-07-23, but profile was on 2019-07-22
water.column.data[which(water.column.data$RM == 318 &
                          water.column.data$date == "2019-07-22"), "date"] <- as.Date("2019-07-23")
# Set the height of the deepest water column sample to an elevation of 1.25, so we can plot it above the corewater profiles
water.column.data <- water.column.data %>%
  full_join(unique.sample.locations) %>%
  filter(corewater == "yes") %>%
  group_by(RM, date) %>%
  filter(elevation_m == min(elevation_m)) %>%
  mutate(height = 1.25) %>%
  select(RM, date, height,
         FTHG, FMHG)



#### Prepare porewater data for sites where we have it ####
porewater.data <- read_xlsx("dataRaw/dataRelease_sed/T3_Shallow sed porewater_HCC.2014-18.xlsx",
                               sheet = "T3_HCC.sed.porewater",
                            skip = 2) %>%
  mutate(RM = `River Mile [Snake R.]`,
         date = as.Date(mdy(`Collection Date (mm/dd/yyyy)`)),
         # Set PW depth to 0.05 m for plotting purposes
         height = -0.05,
         FMHG = `pw.MeHg (ng/L)`,
         FTHG = `pw.THg (ng/L)`) %>%
  group_by(RM, date, height) %>%
  summarise(FTHG = mean(FTHG),
            FMHG = mean(FMHG)) %>%
  ungroup()
  




#### Combine the water column data and porewater data with the corewater profiles ####
all.corewater.data <- rbind(corewater.data,
                            water.column.data,
                            porewater.data)





#### Generate plot function ####
corewater.profile.function <- function(corewater.data.to.use = all.corewater.data,
                                       dates.to.use,
                                       RM.to.use = c(286, 300, 310),
                                       title.to.use) {
  corewater.data.to.use %>%
    filter(date %in% dates.to.use,
           RM %in% RM.to.use) %>%
    ggplot(aes(x = height,
               y = FMHG)) +
    geom_point() +
    # geom_line() +
    facet_wrap(~RM,
               nrow = 1) +
    theme_bw() +
    coord_flip() +
    geom_vline(xintercept = 1.1,
               linetype = "dashed",
               color = cb.translator["vermillion"],
               size = 1) +
    geom_vline(xintercept = 0,
               linetype = "dotted",
               color = cb.translator["reddishpurple"],
               size = 1) +
    xlim(c(-0.1,1.3)) +
    scale_y_continuous(limits = c(0.01, 3.5),
                       trans = 'log10') +
    labs(title = title.to.use,
         y = "Filter-passing MeHg (ng/L)",
         x = "Height above sediment (m)")
  
}



#### Generate all plots
oct.2016 <- corewater.profile.function(dates.to.use = as.Date(c("2016-10-03", "2016-10-05", "2016-10-06")),
                                       title.to.use = "2016 October profiles")
june.2017 <- corewater.profile.function(dates.to.use = as.Date(c("2017-06-05", "2017-06-06", "2017-06-08")),
                                        title.to.use = "2017 June profiles")
sept.2017 <- corewater.profile.function(dates.to.use = as.Date(c("2017-09-25", "2017-09-26", "2017-09-28")),
                                        title.to.use = "2017 September profiles")

june.2018 <- corewater.profile.function(dates.to.use = as.Date(c("2018-06-18", "2018-06-19")),
                                        title.to.use = "2018 June profiles")
sept.2018 <- corewater.profile.function(dates.to.use = as.Date(c("2018-09-24", "2018-09-25", "2018-09-26")),
                                        title.to.use = "2018 September profiles")
july.2019 <- corewater.profile.function(dates.to.use = as.Date(c("2019-07-23", "2019-07-24", "2019-07-25")),
                                        title.to.use = "2019 July profiles")


#### Generate all plots together ####
pdf("results/geochem/MeHg_corewater.pdf",
    height = 14,
    width = 5)
ggarrange(oct.2016, june.2017, sept.2017, june.2018, sept.2018, july.2019,
          ncol = 1)
dev.off()

