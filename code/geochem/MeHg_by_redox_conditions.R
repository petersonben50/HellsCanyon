#### code/geochem/MeHg_by_redox_conditions.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(tidyverse)
source("code/HCC_plotting_needs.R")

#### Read in data ####
geochem.data.adj <- readRDS("dataEdited/geochem/geochem_WC_adjusted_for_redox_classification.rds") %>%
  filter(year(date) != 2015) %>%
  mutate(redox_status = fct_relevel(redox_status, names(renaming.vector)))


renaming.vector <- gsub(", ", "\n", renaming.vector)


#### Generate plot ####
MeHg.redox.boxplot <- geochem.data.adj %>%
  ggplot(aes(x = redox_status,
             y = MeHg_diss_ngL)) +
  geom_boxplot() +
  geom_jitter(aes(shape = as.character(year(date)),
                  col = redox_status),
              width = 0.2) +
  scale_color_manual(values = color.vector) +
  scale_shape_manual(values = shape.vector,
                     labels = renaming.vector,
                     name = "Year") +
  theme_classic() +
  scale_y_continuous(limits = c(0.009, 10),
                     trans = 'log10') +
  scale_x_discrete(labels = renaming.vector) +
  labs(x = "",
       y = "Dissolved MeHg (ng/L)") +
  theme(legend.position = c(0.6, 0.3),
        legend.box.background = element_rect(colour = "black",
                                             size = 1.2),
        legend.key.size = unit(5, units = "mm"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) +
  geom_hline(yintercept = 0.01,
             linetype = 2) +
  annotate(
    "text", label = "Detection limit",
    x = 5, y = 0.012, size = 4, colour = "black"
  )
  # geom_text(aes(x = 5,
  #               y = 0.012,
  #               label = "Detection limit"))
MeHg.redox.boxplot


#### Save out plot ####
pdf("results/geochem/MeHg_by_redoxStatus.pdf",
    height = 5,
    width = 7.5)
MeHg.redox.boxplot
dev.off()


#### Check out averages, SD ####
geochem.data.adj %>%
  group_by(redox_status) %>%
  summarize(MeHg_averages = mean(MeHg_diss_ngL),
            MeHg_max = max(MeHg_diss_ngL),
            MeHg_SD = sd(MeHg_diss_ngL),
            count = n())
