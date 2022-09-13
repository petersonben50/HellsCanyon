#### code/geochem/MeHg_by_redox_conditions.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")

geochem.data.adj <- readRDS("dataEdited/geochem/geochem_WC_adjusted_for_redox_classification.rds")


unique(geochem.data.adj$redox_status)
redox.color.vector <- cb.translator[c("bluishgreen", "reddishpurple", "orange", "skyblue", "blue")]
names(redox.color.vector) <- c("oxic", "suboxic", "no_nitrate_no_sulfide", "no_nitrate_possible_sulfide", "sulfidic")

redox.shape.vector <- c(16, 17, 15, 5, 18)
names(redox.shape.vector) <- names(redox.color.vector)



geochem.data.adj %>%
  # filter(year(date) == 2019) %>%
  # filter(redox_status != "oxic") %>%
  ggplot(aes(x = f_mn_mg_per_l,
             y = MeHg_diss_ngL,
             color = redox_status,
             shape = redox_status)) +
  geom_point() +
  scale_color_manual(values = redox.color.vector,
                     name = "Redox status") +
  scale_shape_manual(values = redox.shape.vector,
                     name = "Redox status") +
  theme_classic() +
  scale_x_continuous(limits = c(0.0001, 2),
                     trans = 'log10') +
  scale_y_continuous(limits = c(0.01, 4),
                     trans = 'log10') +
  labs(x = "Filter-passing Mn (mg/L)",
       y = "Dissolved MeHg (ng/L)") +
  theme(legend.position = c(0.4, 0.88))


