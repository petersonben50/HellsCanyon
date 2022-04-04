#### code/hgcA_analysis/MeHg_hgcA_profiles_SI_figure.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Read in data ####
Hg.data <- read.csv("dataEdited/waterChemistry/Hg_data.csv") %>%
  group_by(date, RM, depth, constituent) %>%
  summarise(concentration = sum(concentration))
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


#### Set up vectors for Hg and hgcA depth plots ####
color.vector <- c(cb.translator["black"],
                  cb.translator["bluishgreen"])
names(color.vector) <- c("hgcA_coverage", "FMHG")
labels.vector <- c("hgcA coverage\n(per 100X SCG)",
                   "Dissolved MeHg\n(ng/L)")
names(labels.vector) <- c("hgcA_coverage", "FMHG")
points.vector <- c(16, 17)
names(points.vector) <- c("hgcA_coverage", "FMHG")


#### Prepare data for plotting ####
MeHg.hgcA.data <- all.data %>%
  filter((year(date) %in% c(2017, 2018) & RM %in% c(286, 300)) |
           (year(date) == 2019 & RM %in% c(300, 310))) %>%
  spread(key = constituent,
         value = concentration) %>%
  filter(!is.na(hgcA_coverage)) %>%
  gather(key = constituent,
         value = abundance,
         -c(1:4)) %>%
  filter(constituent %in% names(color.vector)) %>%
  mutate(date.site = paste(year(date), "-RM", RM,
                           sep = ""))

### Generate depth plots ####
depth.plots <- MeHg.hgcA.data %>%
  ggplot(aes(y = abundance,
             x = depth,
             group = constituent,
             shape = constituent)) +
  geom_point(aes(color = constituent)) +
  geom_line(aes(color = constituent)) +
  scale_colour_manual(values = color.vector,
                      labels = labels.vector) +
  facet_wrap(~date.site, nrow = 3) +
  scale_shape_manual(values = points.vector,
                     labels = labels.vector) +
  # scale_linetype_manual(values = line.vector,
  #                       labels = labels.vector) +
  coord_flip(xlim = c(80, 0)) +
  theme_classic() +
  theme(legend.position = c(0.33, 0.90),
        legend.title = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) +
  xlab("Depth (m)") +
  ylab("Concentration (ng/L)\nGene coverage (per 100X SCG coverage")


pdf("results/manuscript_figures/MeHg_hgcA_profiles_SI_figure.pdf",
    width = 3.2,
    height = 6)
depth.plots
dev.off()
