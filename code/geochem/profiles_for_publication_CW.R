#### code/geochem/profiles_for_publication_CW.R ####
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



#### Set up variables ####
Mn.plotting.factor = 5
Mn.color = "orange"
nitrate.plotting.factor = 5
nitrate.color = "reddishpurple"
MeHg.plotting.factor = 2.5
MeHg.color = "bluishgreen"
min.height = -15
max.height = 110
size.of.points = 1.5
line.width = 2.5



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





corewater.profile.MeHg.and.HgT.function <- function(corewater.data.to.use = corewater.data,
                                                    date.range = c("2016-10-03","2016-10-06"),
                                                    RMs.to.use = c(286, 300, 310),
                                                    title.to.use = "MeHg corewater profiles",
                                                    min.height = -10,
                                                    max.height = 75) {
  corewater.data.to.use <- corewater.data.to.use %>%
    filter(RM %in% RMs.to.use,
           (date >= as.Date(date.range[1]) &
             date <= as.Date(date.range[2])),
           !is.na(FMHG)) %>%
    mutate(height_cm = height * 100)
  
  porewater.data.to.use <- porewater.data %>%
    filter(RM %in% RMs.to.use,
           (date >= as.Date(date.range[1]) &
              date <= as.Date(date.range[2]))) %>%
    mutate(height_cm = -5)
    

  for (RM.of.interest in RMs.to.use) {
    
    corewater.data.site <- corewater.data.to.use %>%
      filter(RM == RM.of.interest) %>%
      arrange(height_cm)
    porewater.data.site <- porewater.data.to.use %>%
      filter(RM == RM.of.interest)
    date.of.sampling.at.RM <- corewater.data.site %>%
      mutate(date = as.character(date)) %>%
      select(date) %>%
      unlist(use.names = FALSE) %>%
      unique()
    
    
    #### Set up plot ####
    plot(x = NULL,
         y = NULL,
         xlab = "",
         ylab = "Elevation (cm)",
         xlim = c(0, 11),
         ylim = c(min.height, max.height),
         xaxt = "n",
         yaxt = "n")
    
    title(main = paste(date.of.sampling.at.RM, ": RM", RM.of.interest, sep = ""),
          cex.main = 0.8)
    # Add depth axis
    axis(2,
         at = seq(0, 100, by = 50),
         labels = seq(0, 100, by = 50),
         las = 2,
         cex.axis = 0.85)
    
    
    #### Add Mn ####
    lines(x = corewater.data.site$f_mn_mg_per_l*Mn.plotting.factor,
          y = corewater.data.site$height_cm,
          col = cb.translator[Mn.color],
          lwd = line.width)
    points(x = corewater.data.site$f_mn_mg_per_l*Mn.plotting.factor,
           y = corewater.data.site$height_cm,
           bg = cb.translator[Mn.color],
           pch = 24,
           cex = size.of.points)
    # Add axis for Mn
    axis(1,
         at = seq(0, 10, by = 2.5),
         labels = seq(0, 10, by = 2.5)/Mn.plotting.factor,
         cex.axis = 0.85)
    # Add title for Mn
    title(xlab = "Mn (mg/L)",
          line = 1.2,
          cex.lab = 0.9)
    
    #### Add MeHg ####
    lines(x = corewater.data.site$FMHG*MeHg.plotting.factor,
          y = corewater.data.site$height_cm,
          col = cb.translator[MeHg.color],
          lwd = line.width)
    points(x = corewater.data.site$FMHG*MeHg.plotting.factor,
           y = corewater.data.site$height_cm,
           bg = cb.translator[MeHg.color],
           pch = 24,
           cex = size.of.points)

    #### Add PW MeHg and Mn ####
    points(x = porewater.data.site$FMHG*MeHg.plotting.factor,
           y = porewater.data.site$height_cm,
           col = cb.translator[MeHg.color],
           pch = 5,
           cex = size.of.points)
    points(x = porewater.data.site$f_mn_mg_per_l*Mn.plotting.factor,
           y = porewater.data.site$height_cm,
           col = cb.translator[Mn.color],
           pch = 5,
           cex = size.of.points)
    
    #### Add sediment-water interface line ####
    abline(h = 0)
    
  }
}


pdf("results/geochem/profiles_for_figure/CW_profiles.pdf",
    width = 7.5,
    height = 2)

par(mfrow = c(2,6),
    mar = c(2, 2.5, 1, 0.5),
    mgp=c(1.5,0.2,0),
    tck=-0.008)


corewater.profile.MeHg.and.HgT.function(date.range = c("2016-10-03","2016-10-06"),
                                        max.height = 110)
corewater.profile.MeHg.and.HgT.function(date.range = c("2017-09-25", "2017-09-28"),
                                        max.height = 110)
corewater.profile.MeHg.and.HgT.function(date.range = c("2018-09-24", "2018-09-26"),
                                        max.height = 110)
corewater.profile.MeHg.and.HgT.function(date.range = c("2019-07-22", "2019-07-25"),
                                        max.height = 110)
dev.off()
