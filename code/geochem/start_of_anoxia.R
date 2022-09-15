
setwd("~/Documents/research/HellsCanyon/")
library(fields)
library(lubridate)
library(plotly)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")

#### Read in data ####
seabird.data <- read_xlsx("dataRaw/dataRelease_chem/v2/HCC - V2 Data Release_Master File_Oct 2019 to present_V2_9_forModelingGroup.xlsx",
                         sheet = "Table_6_Brownlee_Profiles") %>%
  rename(depth = measurement_depth_m,
         elevation_m = brownlee_reservoir_measurement_elevation_m) %>%
  mutate(RM = as.character(snake_river_mile),
         date = as.Date(measurement_collection_date_mm_dd_yy_h_mm)) %>%
  select(RM, date, depth, elevation_m, water_temp_deg_c, diss_oxy_mg_per_l,
         spec_cond_ms_per_cm, ph, turb_ntu)



#### Generate matrix with min DO value for each profile (by date and RM) ####
min.do.data <- seabird.data %>%
  group_by(date, RM) %>%
  summarize(DO = min(diss_oxy_mg_per_l)) %>%
  ungroup() %>%
  arrange(date, RM) %>%
  mutate(DO = DO + 0.001)

profiling.dates.to.include <- min.do.data %>%
  group_by(date) %>%
  summarise(profiles_taken = n()) %>%
  filter(profiles_taken >= 5) %>%
  select(date) %>%
  mutate(date = as.character(date)) %>%
  unlist(use.names = FALSE)
 

#### Function to generate contour plots of DO ####

plot.anoxia.by.year <- function(year.of.interest) {
  min.do.data.year.of.interest <- min.do.data %>%
    filter(year(date) == year.of.interest) %>%
    spread(key = RM,
           value = DO)
    


  min.do.data.matrix <- min.do.data.year.of.interest %>%
    select(-date) %>%
    as.matrix() %>%
    t() %>%
    image.smooth()
  min.do.data.matrix$y <- min.do.data.year.of.interest$date
  min.do.data.matrix$x <- as.numeric(colnames(min.do.data.year.of.interest)[-1])

  filled.contour(x = min.do.data.matrix$x,
                 y = min.do.data.matrix$y,
                 z = min.do.data.matrix$z,
                 ylim = c(as.Date(paste(year.of.interest, "-12-01",
                                        sep = "")),
                          as.Date(paste(year.of.interest, "-05-01",
                                        sep = ""))),
                 xlim = c(286, 330),
                 zlim = c(0, 10),
                 nlevels = 40,
                 color.palette = colorRampPalette(c(cb.translator["black"], cb.translator["blue"], cb.translator["bluishgreen"], cb.translator["orange"]),
                                                  bias = 1, space = "rgb"))

}



#### Save out PDF of plots ####
pdf("results/geochem/start_of_anoxia.pdf",
    width = 6,
    height = 6)
plot.anoxia.by.year(2015)
plot.anoxia.by.year(2016)
plot.anoxia.by.year(2017)
plot.anoxia.by.year(2018)
dev.off()


