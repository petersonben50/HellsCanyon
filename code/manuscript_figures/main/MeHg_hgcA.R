#### code/manuscript_figures/MeHg_hgcA.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Read in data ####
Hg.data <- read.csv("dataEdited/waterChemistry/Hg_data.csv")
hgcA.data <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv")


#### Summarize hgcA data at each depth ####
hgcA.data <- hgcA.data %>%
  filter(rep == TRUE) %>%
  filter(!(RM %in% c(305, 314, 318)),
         !(RM == 286 & year(date) == 2019)) %>%
  group_by(date, RM, depth, redoxClassification) %>%
  summarise(hgcA_coverage = sum(coverage)) %>%
  ungroup()


#### Combine data ####
all.data <- full_join(Hg.data,
                      hgcA.data)
rm(Hg.data,
   hgcA.data)


#### Make Hg and hgcA depth plots ####
color.vector <- c(cb.translator["black"],
                  cb.translator["bluishgreen"])
names(color.vector) <- c("hgcA_coverage", "FMHG")
labels.vector <- c("hgcA coverage\n(per 100X SCG)",
                  "Dissolved MeHg\n(ng/L)")
names(labels.vector) <- c("hgcA_coverage", "FMHG")
points.vector <- c(16, 17)
names(points.vector) <- c("hgcA_coverage", "FMHG")

depth.plots <- all.data %>%
  filter(year(date) %in% c(2017, 2018),
         RM %in% c(286, 300)) %>%
  gather(key = constituent,
         value = abundance,
         -c(1:3, 8)) %>%
  filter(!is.na(abundance),
         constituent %in% names(color.vector)) %>%
  mutate(date.site = paste(year(date), "-RM", RM,
                           sep = "")) %>%
  ggplot(aes(y = abundance,
             x = depth,
             group = constituent,
             shape = constituent)) +
  geom_point(aes(color = constituent)) +
  geom_line(aes(color = constituent)) +
  scale_colour_manual(values = color.vector,
                      labels = labels.vector) +
  facet_wrap(~date.site, nrow = 2) +
  scale_shape_manual(values = points.vector,
                     labels = labels.vector) +
  # scale_linetype_manual(values = line.vector,
  #                       labels = labels.vector) +
  coord_flip(xlim = c(80, 0)) +
  theme_classic() +
  theme(legend.position = c(0.83, 0.90),
        legend.title = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) +
  xlab("Depth (m)") +
  ylab("Concentration (ng/L)\nGene coverage (per 100X SCG coverage")


#### Make MeHg vs hgcA scatterplot ####
redox.color.vector <- c(cb.translator["bluishgreen"],
                        cb.translator["orange"],
                        cb.translator["blue"])
names(redox.color.vector) <- c("oxic", "suboxic", "euxinic")
redox.year.vector <- c(16, 17, 2)
names(redox.year.vector) <- c("2017", "2018", "2019")
hgcA.MeHg.scatterplot <- all.data %>%
  filter(!is.na(hgcA_coverage)) %>%
  mutate(redoxClassification = fct_relevel(redoxClassification,
                                           names(redox.color.vector))) %>%
  ggplot(aes(x = log(hgcA_coverage, 10),
             y = log(FMHG, 10),
             group = redoxClassification,
             shape = as.character(year(date)))) +
  geom_point(aes(color = redoxClassification)) +
  scale_color_manual(values = redox.color.vector,
                     name = "Redox status") +
  scale_shape_manual(values = redox.year.vector,
                     name = "Year") +
  xlab("Log of hgcA abundance") +
  ylab("Log of dissolved MeHg (ng/L)") +
  theme_classic() +
  theme(legend.position = c(0.22, 0.75),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.key = element_rect(fill = "transparent", colour = "black"))



pdf("results/manuscript_figures/MeHg_hgcA.pdf",
    width = 8,
    height = 5)
depth.plots + hgcA.MeHg.scatterplot
dev.off()
