#### code/binning/depth_plots.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
cb.translator <- c(cb.translator, "grey50")
names(cb.translator)[length(cb.translator)] <- "grey"


#### Read in metadata ####
MG.metadata <- read_xlsx("metadata/metagenome_metadata.xlsx")


#### Read in normalization table ####
MG.metadata <- read_xlsx("metadata/metagenome_metadata.xlsx")



#### Read in metabolism data ####
metabolic.data <- read_xlsx("dataEdited/binning/metabolism/metabolic_summary.xlsx") %>%
  select(HMS, binID, metabolic_assignment)
HMS.colors <- read_xlsx("dataEdited/binning/metabolism/metabolic_summary.xlsx",
                        sheet = "HMS_colors") %>%
  select(HMS, colorToUse)
HMS.colors.vector <- cb.translator[HMS.colors$colorToUse]
names(HMS.colors.vector) <- HMS.colors$HMS


#### Read in depth data ####
depth.2017 <- read.table("dataEdited/binning/coverageAnvio/coverage_goodBins_2017.txt",
                         header = TRUE) %>%
  rename(binID = bins) %>%
  gather(key = metagenomeID,
         value = coverage,
         -1) %>%
  mutate(metagenomeID = metagenomeID %>% strsplit("_") %>% sapply("[", 1)) %>%
  left_join(MG.metadata %>% select(metagenomeID, RM, depth)) %>%
  left_join(metabolic.data)

depth.2017 %>%
  ggplot(aes(y = coverage,
             x = depth,
             group = binID)) +
  geom_point(aes(color = HMS)) +
  geom_line(aes(color = HMS)) +
  scale_color_manual(values = HMS.colors.vector) +
  facet_wrap(~metabolic_assignment + RM, nrow = 4) +
  coord_flip() +
  scale_x_reverse(c(0, 80)) +
  theme_bw()



#### Read in depth data ####
depth.2018 <- read.table("dataEdited/binning/coverageAnvio/coverage_goodBins_2018.txt",
                         header = TRUE) %>%
  rename(binID = bins) %>%
  gather(key = metagenomeID,
         value = coverage,
         -1) %>%
  mutate(metagenomeID = metagenomeID %>% strsplit("_") %>% sapply("[", 1)) %>%
  left_join(MG.metadata %>% select(metagenomeID, RM, depth)) %>%
  left_join(metabolic.data) %>%
  arrange(depth)

depth.2018 %>%
  ggplot(aes(y = coverage,
             x = depth,
             group = binID)) +
  geom_point(aes(color = HMS)) +
  geom_line(aes(color = HMS)) +
  scale_color_manual(values = HMS.colors.vector) +
  facet_wrap(~metabolic_assignment + RM, nrow = 3) +
  coord_flip() +
  scale_x_reverse(c(0, 80)) +
  theme_bw()



#### Read in depth data ####
depth.2019 <- read.table("dataEdited/binning/coverageAnvio/coverage_goodBins_2019.txt",
                         header = TRUE) %>%
  rename(binID = bins) %>%
  gather(key = metagenomeID,
         value = coverage,
         -1) %>%
  mutate(metagenomeID = metagenomeID %>% strsplit("_") %>% sapply("[", 1)) %>%
  left_join(MG.metadata %>% select(metagenomeID, RM, depth)) %>%
  filter(RM %in% c(300, 310)) %>%
  left_join(metabolic.data) %>%
  arrange(depth)

depth.2019 %>%
  ggplot(aes(y = coverage,
             x = depth,
             group = binID)) +
  geom_point(aes(color = HMS)) +
  geom_line(aes(color = HMS)) +
  scale_color_manual(values = HMS.colors.vector) +
  facet_wrap(~metabolic_assignment + RM, nrow = 2) +
  coord_flip() +
  scale_x_reverse(c(0, 80)) +
  theme_bw()
