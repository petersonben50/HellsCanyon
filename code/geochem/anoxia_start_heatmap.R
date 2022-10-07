#### code/geochem/anoxia_start_heatmap.R ####
# Benjamin D. Peterson

#### Set the table ####
setwd("~/Documents/research/HellsCanyon/")
library(fields)
library(lubridate)
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


#### Cap maximum DO (to highlight changes at lower concentrations) ####
min.do.data$DO[which(min.do.data$DO > 10)] <- 10


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

  # Calculate date and location of anoxia onset
  date.of.anoxia <- min.do.data %>%
    filter(year(date) == year.of.interest) %>%
    arrange(date) %>%
    filter(DO < 0.5) %>%
    filter(date == min(date))
  
  # Generate plot
  filled.contour(x = min.do.data.matrix$x,
                 y = min.do.data.matrix$y,
                 z = min.do.data.matrix$z,
                 plot.axes={
                   axis(1);
                   axis(2,
                        at = floor_date(seq(from = as.Date(paste(year.of.interest, "-04-01",
                                                                 sep = "")),
                                            to = as.Date(paste(year.of.interest, "-11-01",
                                                               sep = "")),
                                            by = 1),
                                        'month') %>%
                          unique(),
                        labels = month(seq(from = as.Date(paste(year.of.interest, "-04-01",
                                                                     sep = "")),
                                                to = as.Date(paste(year.of.interest, "-11-01",
                                                                   sep = "")),
                                                by = 1),
                                       label = TRUE,
                                       abbr = TRUE) %>%
                          unique());
                   points(x = date.of.anoxia$RM,
                          y = date.of.anoxia$date,
                          col = cb.translator["reddishpurple"],
                          pch = 16)},
                 ylim = c(as.Date(paste(year.of.interest, "-11-01",
                                        sep = "")),
                          as.Date(paste(year.of.interest, "-04-01",
                                        sep = ""))),
                 xlim = c(286, 325),
                 zlim = c(0, 10),
                 plot.title = title(main = year.of.interest),
                 nlevels = 40,
                 color.palette = colorRampPalette(c(cb.translator["black"], cb.translator["blue"], cb.translator["bluishgreen"], cb.translator["orange"]),
                                                  bias = 1, space = "rgb"))

}


#### Save out PDF of plots ####
par(mar = c(3, 3, 1, 1),
    mgp=c(1.5,0.4,0),
    tck=-0.008)
pdf("results/geochem/start_of_anoxia.pdf",
    width = 6,
    height = 5)
plot.anoxia.by.year(2015)
plot.anoxia.by.year(2016)
plot.anoxia.by.year(2017)
plot.anoxia.by.year(2018)
plot.anoxia.by.year(2019)
dev.off()


