#### code/hgcA_analyses/hgcA_vs_MeHg.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon")
library(lubridate)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Read in hgcA data ####
hgcA.data <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv") %>%
  filter(rep) %>%
  group_by(RM, depth, date, metagenomeID) %>%
  summarise(coverage = sum(coverage))


#### Read in MeHg data ####
Hg.data.2017.2018 <- read.csv("dataEdited/waterChemistry/Hg_2015_2018.csv",
                              stringsAsFactors = FALSE) %>%
  filter(year(date) %in% c("2017", "2018")) %>%
  filter(month(date) == "9") %>%
  filter(day(date) >= 23) %>%
  select(RM, depth, date, FTHG,
         FMHG, PTHG, PMHG)
Hg.data.2019 <- read.csv("dataEdited/waterChemistry/Hg_WC_2019_intensive.csv",
                         stringsAsFactors = FALSE) %>%
  select(RM, depth, date, FTHG,
         FMHG, PTHG, PMHG)
Hg.data <- rbind(Hg.data.2017.2018, Hg.data.2019) %>%
  gather(key = constituent,
         value = concentration,
         4:7) %>%
  group_by(RM, depth, date, constituent) %>%
  summarise(concentration = mean(concentration)) %>%
  ungroup() %>%
  spread(key = constituent,
         value = concentration)
rm(Hg.data.2017.2018,
   Hg.data.2019)


#### Read in sulfide data ####
geochem.data <- read.csv("dataEdited/waterChemistry/geochem_WC_2015_2018.csv",
                         stringsAsFactors = FALSE) %>%
  filter(year(date) %in% c("2017", "2018")) %>%
  filter(month(date) == "9") %>%
  select(RM, depth, date, sulfide_mg.L)
# Replace non-detects with zeroes in sulfide column
geochem.data$sulfide_mg.L[geochem.data$sulfide_mg.L == "nd"] <- 0
geochem.data$sulfide_mg.L <- as.numeric(geochem.data$sulfide_mg.L)



#### Combine data ####
all.data <- hgcA.data %>%
  left_join(Hg.data) %>%
  left_join(geochem.data) %>%
  mutate(year.sampling = as.character(year(date))) %>%
  mutate(point_shape = "16")
all.data[all.data$sulfide_mg.L > 0, "point_shape"] <- "18"

#### Set up color vector ####
color.vector.year <- c(cb.translator["skyblue"],
                       cb.translator["bluishgreen"],
                       cb.translator["vermillion"])
names(color.vector.year) <- c("2017", "2018", "2019")


#### Plot dissolved MeHg against hgcA coverage ####
pdf("results/hgcA_analysis/MeHg_vs_hgcA.pdf",
    height = 4,
    width = 4)
all.data %>%
  ggplot(aes(x = log(coverage),
             y = log(FMHG),
             fill = year.sampling)) +
  geom_point(aes(color = year.sampling)) +
  scale_color_manual(values = color.vector.year) +
  labs(x = "log(hgcA abundance)",
       y = "log(Dissolved MeHg (ng/L))") +
  theme_classic() +
  theme(legend.position = c(0.2, 0.8),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, "mm"),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
dev.off()


#### Combine with sulfide ####
all.data.sulfide <- all.data %>%
  left_join(geochem.data) %>%
  filter(year.sampling != 2019) %>%
  mutate(point_shape = "sulfide-n.d.")
all.data.sulfide[all.data.sulfide$sulfide_mg.L > 0, "point_shape"] <- "sulfide-detected"
shape.vector <- c(16, 18)
names(shape.vector) <- c("sulfide-n.d.", "sulfide-detected")


#### Plot dissolved MeHg against hgcA coverage without 2019 data ####
pdf("results/hgcA_analysis/MeHg_vs_hgcA_2017_2018.pdf",
    height = 4,
    width = 4)
all.data.sulfide %>%
  # filter(year.sampling != "2019") %>%
  ggplot(aes(x = log(coverage),
             y = log(FMHG),
             fill = year.sampling)) +
  geom_point(aes(color = year.sampling,
                 shape = point_shape)) +
  scale_color_manual(values = color.vector.year) +
  scale_shape_manual(values = shape.vector) +
  labs(x = "log(hgcA abundance)",
       y = "log(Dissolved MeHg (ng/L))") +
  theme_classic() +
  theme(legend.position = c(0.2, 0.8),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.spacing.y = unit(0, "mm"))
dev.off()
