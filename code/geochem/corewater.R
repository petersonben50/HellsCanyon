#### code/geochem/corewater.R ####
# Written by Benjamin D. Peterson


#### Get ready. Get set ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(ggpubr)
library(grid)
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
         f_mn_mg_per_l = (as.numeric(f_mn_mcg_per_l) / 1000),
         FTHG = as.numeric(f_thg_ng_per_l),
         FMHG = as.numeric(f_mehg_ng_per_l)) %>%
  mutate(height = (upper_height + lower_height) / 2) %>%
  select(RM, date, height, FTHG, FMHG,
         f_inorganic_sulfide_mg_per_l,
         f_mn_mg_per_l) %>%
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
         f_mn_mg_per_l = (as.numeric(f_mn_mcg_per_l) / 1000),
         FTHG = as.numeric(f_thg_ng_per_l),
         FMHG = as.numeric(f_mehg_ng_per_l)) %>%
  filter(RM != "--",
         depth != "--") %>%
  select(RM, date, depth, elevation_m, FTHG, FMHG,
         f_inorganic_sulfide_mg_per_l,
         f_mn_mg_per_l) %>%
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
         FTHG, FMHG,
         f_inorganic_sulfide_mg_per_l,
         f_mn_mg_per_l)



#### Prepare porewater data for sites where we have it ####
porewater.data <- read_xlsx("dataRaw/dataRelease_sed/T3_Shallow sed porewater_HCC.2014-18.xlsx",
                               sheet = "T3_HCC.sed.porewater",
                            skip = 2) %>%
  mutate(RM = `River Mile [Snake R.]`,
         date = as.Date(mdy(`Collection Date (mm/dd/yyyy)`)),
         # Set PW depth to 0.05 m for plotting purposes
         height = -0.05,
         FMHG = `pw.MeHg (ng/L)`,
         FTHG = `pw.THg (ng/L)`,
         f_inorganic_sulfide_mg_per_l = as.numeric(`pw.S(-II) (mg/L)`),
         f_mn_mg_per_l = as.numeric(`pw.Mn (µg/L)`)/1000) %>%
  group_by(RM, date, height) %>%
  summarise(FTHG = mean(FTHG),
            FMHG = mean(FMHG),
            f_inorganic_sulfide_mg_per_l = mean(f_inorganic_sulfide_mg_per_l),
            f_mn_mg_per_l = mean(f_mn_mg_per_l)) %>%
  ungroup()
  




#### Combine the water column data and porewater data with the corewater profiles ####
all.corewater.data <- rbind(corewater.data,
                            water.column.data,
                            porewater.data)



#### Add early vs. late sampling time info ####
all.corewater.data <- all.corewater.data %>%
  mutate(fall.sample = (month(date) >8))
all.corewater.data[which(all.corewater.data$fall.sample), "sampling.period"] <- "Late Stratification"
all.corewater.data[which(!all.corewater.data$fall.sample), "sampling.period"] <- "Early Stratification"
all.corewater.data <- all.corewater.data %>%
  select(RM, date, height, FTHG, FMHG,
         sampling.period,
         f_inorganic_sulfide_mg_per_l,
         f_mn_mg_per_l)
  
  
#### Generate plot function for both MeHg and HgT ####
# corewater.profile.MeHg.and.HgT.function <- function(corewater.data.to.use = all.corewater.data,
#                                                     dates.to.use,
#                                                     RM.to.use = c(286, 300, 310),
#                                                     title.to.use) {
#   corewater.data.to.use %>%
#     filter(date %in% dates.to.use,
#            RM %in% RM.to.use) %>%
#     ggplot() +
#     geom_point(aes(x = height,
#                    y = FMHG)) +
#     geom_point(aes(x = height,
#                    y = FTHG), 
#                col = cb.translator["vermillion"]) +
#     facet_wrap(~RM,
#                nrow = 1) +
#     theme_bw() +
#     coord_flip() +
#     geom_vline(xintercept = 1.1,
#                linetype = "dashed",
#                color = cb.translator["vermillion"],
#                size = 1) +
#     geom_vline(xintercept = 0,
#                linetype = "dotted",
#                color = cb.translator["reddishpurple"],
#                size = 1) +
#     xlim(c(-0.1,1.3)) +
#     scale_y_continuous(limits = c(0.01, 45),
#                        trans = 'log10') +
#     labs(title = title.to.use,
#          y = "Filter-passing MeHg (ng/L)",
#          x = "Height above sediment (m)")
#   
# }
# 
# 
# 
# #### Generate all plots
# oct.2016 <- corewater.profile.MeHg.and.HgT.function(dates.to.use = as.Date(c("2016-10-03", "2016-10-05", "2016-10-06")),
#                                        title.to.use = "2016 October profiles")
# june.2017 <- corewater.profile.MeHg.and.HgT.function(dates.to.use = as.Date(c("2017-06-05", "2017-06-06", "2017-06-08")),
#                                         title.to.use = "2017 June profiles")
# sept.2017 <- corewater.profile.MeHg.and.HgT.function(dates.to.use = as.Date(c("2017-09-25", "2017-09-26", "2017-09-28")),
#                                         title.to.use = "2017 September profiles")
# 
# june.2018 <- corewater.profile.MeHg.and.HgT.function(dates.to.use = as.Date(c("2018-06-18", "2018-06-19")),
#                                         title.to.use = "2018 June profiles")
# sept.2018 <- corewater.profile.MeHg.and.HgT.function(dates.to.use = as.Date(c("2018-09-24", "2018-09-25", "2018-09-26")),
#                                         title.to.use = "2018 September profiles")
# july.2019 <- corewater.profile.MeHg.and.HgT.function(dates.to.use = as.Date(c("2019-07-23", "2019-07-24", "2019-07-25")),
#                                         title.to.use = "2019 July profiles")
# 
# 
# #### Generate all MeHg and HgT plots together ####
# ggarrange(oct.2016, june.2017, sept.2017, june.2018, sept.2018, july.2019,
#           ncol = 1)



#### Generate plot function for multiple years of MeHg data on one plot ####
corewater.profile.MeHg.and.HgT.function <- function(corewater.data.to.use = all.corewater.data,
                                                    RM.to.use = c(286, 300, 310, 318),
                                                    title.to.use = "MeHg corewater profiles") {
  # Text annotations
  water.column.annot <- grobTree(textGrob("Deepest water\nprofile sample",
                                          x = 0.1, y = 0.93, hjust = 0,
                                          gp = gpar(col = "black",
                                                    fontsize = 8,
                                                    fontface = "italic")))
  pore.water.annot <- grobTree(textGrob("Porewater\nsample",
                                        x = 0.1, y = 0.08, hjust = 0,
                                        gp = gpar(col = "black",
                                                  fontsize = 8,
                                                  fontface = "italic")))
  corewater.data.to.use %>%
    filter(RM %in% RM.to.use,
           !is.na(FMHG)) %>%
    ggplot() +
    geom_point(aes(x = height,
                   y = FMHG,
                   col = paste(year(date), month(date, label = TRUE, abbr = TRUE),
                               sep = '-'))) +
    geom_line(aes(x = height,
                  y = FMHG,
                  col = paste(year(date), month(date, label = TRUE, abbr = TRUE),
                              sep = '-'))) + 
    scale_color_manual(values = color.vector.to.use,
                       name = "Sampling Time") +
    facet_wrap(~sampling.period + RM,
               nrow = 2) +
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
    xlim(c(-0.2,1.5)) +
    scale_y_continuous(limits = c(0.01, 10),
                       trans = 'log10') +
    labs(title = title.to.use,
         y = "Filter-passing MeHg (ng/L)",
         x = "Height above sediment (m)") +
    annotation_custom(water.column.annot) +
    annotation_custom(pore.water.annot)
  
}


#### Set up color vector ####
time.points.to.use <- all.corewater.data %>%
  filter(RM %in% c(286, 300, 310, 318),
         !is.na(FMHG)) %>%
  mutate(times.to.use = paste(year(date), month(date, label = TRUE, abbr = TRUE),
                              sep = '-')) %>%
  select(times.to.use) %>%
  unlist(use.names = FALSE) %>%
  unique()
color.vector.to.use = cb.translator[c("bluishgreen", "vermillion", "blue", "orange", "skyblue", "yellow")]
names(color.vector.to.use) <- time.points.to.use


#### Save out pdf of image ####
pdf("results/geochem/MeHg_corewater.pdf",
    height = 7,
    width = 7.2)

corewater.profile.MeHg.and.HgT.function()
dev.off()

