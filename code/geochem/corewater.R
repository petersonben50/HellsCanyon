#### code/geochem/corewater.R ####
# Written by Benjamin D. Peterson


#### Get ready. Get set ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(ggpubr)
library(lubridate)
library(readxl)
library(tidyverse)


#### Read in data ####
geochem.data <- read_xlsx("dataRaw/dataRelease_chem/v2/HCC - V2 Data Release_Master File_Oct 2019 to present_V2_9_forModelingGroup.xlsx",
                          sheet = "Table_3_Water_V2") %>%
  filter(min_sample_height_above_sediment_m != "--") %>%
  mutate(RM = as.character(snake_river_mile),
         date = as.Date(sample_collection_date_mm_dd_yy_h_mm),
         upper_height = as.numeric(max_sample_height_above_sediment_m),
         lower_height = as.numeric(min_sample_height_above_sediment_m),
         FTHG = as.numeric(f_thg_ng_per_l),
         FMHG = as.numeric(f_mehg_ng_per_l)) %>%
  mutate(height = (upper_height + lower_height) / 2) %>%
  select(RM, date, height,
         FTHG, FMHG)



#### Plot 2016 data ####
corewater.2016 <- geochem.data %>%
  filter(year(date) == 2016,
         RM %in% c(286, 300, 310)) %>%
  ggplot(aes(x = height,
             y = FMHG)) +
  geom_point() +
  facet_wrap(~month(date, label = T) + RM,
             nrow = 1) +
  theme_bw() +
  coord_flip() +
  xlim(c(0,1)) +
  scale_y_continuous(limits = c(0.01, 3),
                     trans = 'log10') +
  labs(title = "2016 corewater profiles",
       x = "Filter-passing MeHg (ng/L)",
       y = "Height above sediment (m)")



#### Plot 2017 data ####
corewater.2017 <- geochem.data %>%
  filter(year(date) == 2017,
         RM %in% c(286, 300, 310)) %>%
  ggplot(aes(x = height,
             y = FMHG)) +
  geom_point() +
  facet_wrap(~month(date, label = T) + RM,
             nrow = 2) +
  theme_bw() +
  coord_flip() +
  xlim(c(0,1)) +
  scale_y_continuous(limits = c(0.01, 3),
                     trans = 'log10') +
  labs(title = "2017 corewater profiles",
       x = "Filter-passing MeHg (ng/L)",
       y = "Height above sediment (m)")


#### Plot 2018 data ####
corewater.2018 <- geochem.data %>%
  filter(year(date) == 2018,
         RM %in% c(286, 300, 310)) %>%
  ggplot(aes(x = height,
             y = FMHG)) +
  geom_point() +
  facet_wrap(~month(date, label = T) + RM,
             nrow = 2) +
  theme_bw() +
  coord_flip() +
  xlim(c(0,1)) +
  scale_y_continuous(limits = c(0.01, 3),
                     trans = 'log10') +
  labs(title = "2018 corewater profiles",
       x = "Filter-passing MeHg (ng/L)",
       y = "Height above sediment (m)")


#### Plot 2019 data ####
corewater.2019 <- geochem.data %>%
  filter(year(date) == 2019) %>%
  ggplot(aes(x = height,
             y = FMHG)) +
  geom_point() +
  facet_wrap(~month(date, label = T) + RM,
             nrow = 2) +
  theme_bw() +
  coord_flip() +
  xlim(c(0,1)) +
  scale_y_continuous(limits = c(0.01, 3),
                     trans = 'log10') +
  labs(title = "2019 corewater profiles",
       x = "Filter-passing MeHg (ng/L)",
       y = "Height above sediment (m)")




#### Generate all plots together ####
pdf("results/geochem/MeHg_corewater.pdf",
    height = 14,
    width = 5)
ggarrange(corewater.2016, corewater.2017, corewater.2018, corewater.2019,
          ncol = 1)
dev.off()
