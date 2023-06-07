#### code/hgcA_analysis/hgcA_abundance.R ####
# Written by Benjamin D. Peterson


#### Get ready. Get set ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(lubridate)
library(patchwork)
library(readxl)
library(tidyverse)
source("code/HCC_plotting_needs.R")


#### Read in data ####
hgcA.data <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv") %>%
  filter(rep) %>%
  filter(((RM %in% c(286, 300)) & (year(date) %in% c(2017, 2018))) |
           ((RM %in% c(300, 310)) & (year(date) == 2019))) %>%
  group_by(metagenomeID, date, depth, redoxClassification) %>%
  summarise(coverage = sum(coverage)) %>%
  ungroup() %>%
  mutate(redoxClassification = fct_relevel(redoxClassification,
                                           names(color.vector)[c(1:3, 5)]))
hgcA.data[hgcA.data$coverage < 0.001, "coverage"] <- 0.001

#### Fix up renaming vector ####
renaming.vector <- gsub(", ", "\n", renaming.vector)




#### Check out averages, max, etc for each redox classification ####
hgcA.data %>%
  group_by(redoxClassification) %>%
  summarise(mean_coverage = mean(coverage),
            max_coverage = max(coverage),
            sd_coverage = sd(coverage),
            n_coverage = n())


#### Generate plot ####
hgcA.plot <- hgcA.data %>%
  ggplot(aes(x = redoxClassification,
             y = coverage)) +
  geom_jitter(aes(shape = as.character(year(date)),
                  col = redoxClassification),
              width = 0.2,
              size = 2) +
  scale_color_manual(values = color.vector) +
  scale_shape_manual(values = shape.vector,
                     labels = renaming.vector,
                     name = "Year") +
  theme_classic() +
  scale_y_continuous(limits = c(0.0009, 10),
                     trans = 'log10') +
  scale_x_discrete(labels = renaming.vector) +
  labs(x = "",
       y = "hgcA abundance (%)") +
  theme(legend.position = "bottomright",
        legend.box.background = element_rect(colour = "black",
                                             size = 1.2),
        legend.key.size = unit(5, units = "mm"),
        axis.text.x = element_text(colour = "black",
                                   angle = 45,
                                   vjust = 0.55),
        axis.text.y = element_text(colour = "black"))


#### Save plot ####
pdf("results/hgcA_analysis/abundance_hgcA_by_redox.pdf",
    width = 4,
    height = 4.5)
hgcA.plot
dev.off()
