#### code/geochem/MeHg_by_DOC_SUVA.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(ggpubr)
library(lubridate)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")



geochem.data.adj <- readRDS("dataEdited/geochem/geochem_WC_adjusted_for_redox_classification.rds")




#### Generate needed vectors ####
unique(geochem.data.adj$redox_status)
redox.color.vector <- cb.translator[c("bluishgreen", "reddishpurple", "orange", "black", "blue")]
names(redox.color.vector) <- c("oxic", "suboxic", "no_nitrate_no_sulfide", "no_nitrate_possible_sulfide", "sulfidic")

redox.year.vector <- c(5, 18, 16, 17, 15)
names(redox.year.vector) <- c(2015, 2016, 2017, 2018, 2019)

renaming.vector <- c("Oxygen detected", "No oxygen, nitrate detected", "No nitrate, no sulfide",
                     "No nitrate, sulfide not measured", "Sulfide detected")
names(renaming.vector) <- names(redox.color.vector)


#### Plot DOC by SUVA with redox for overview ####
geochem.data.adj %>%
  filter(!is.na(doc_boulder_mgc_per_l)) %>%
  ggplot(aes(x = doc_boulder_mgc_per_l,
             y = suva_254nm_l_per_mgc_per_m,
             col = redox_status)) +
  geom_point(aes(color = redox_status,
                 shape = as.character(year(date)))) +
  scale_color_manual(values = redox.color.vector,
                     labels = renaming.vector,
                     name = "Redox status") +
  scale_shape_manual(values = redox.year.vector,
                     labels = renaming.vector,
                     name = "Year") +
  theme_classic()


#### Plot MeHg by SUVA with redox for overview ####
MeHg.vs.SUVA <- geochem.data.adj %>%
  filter(!is.na(suva_254nm_l_per_mgc_per_m)) %>%
  ggplot(aes(x = suva_254nm_l_per_mgc_per_m,
             y = MeHg_diss_ngL,
             col = redox_status)) +
  geom_point(aes(color = redox_status,
                 shape = as.character(year(date)))) +
  scale_color_manual(values = redox.color.vector,
                     labels = renaming.vector,
                     name = "Redox status") +
  scale_shape_manual(values = redox.year.vector,
                     labels = renaming.vector,
                     name = "Year") +
  theme_classic() +
  labs(x = "DOC SUVA254",
       y = "MeHg (ng/L)") +
  theme(legend.position = "none")


#### Plot MeHg by DOC with redox for overview ####
MeHg.vs.DOC <- geochem.data.adj %>%
  filter(!is.na(doc_boulder_mgc_per_l)) %>%
  ggplot(aes(x = doc_boulder_mgc_per_l,
             y = MeHg_diss_ngL,
             col = redox_status)) +
  geom_point(aes(color = redox_status,
                 shape = as.character(year(date)))) +
  scale_color_manual(values = redox.color.vector,
                     labels = renaming.vector,
                     name = "Redox status") +
  scale_shape_manual(values = redox.year.vector,
                     labels = renaming.vector,
                     name = "Year") +
  theme_classic() +
  labs(x = "DOC (mg/L)",
       y = "MeHg (ng/L)") +
  theme(legend.position = c(0.3, 0.85))


#### Save out plot ####
pdf("results/geochem/MeHg_DOC.pdf",
    width = 7.2,
    height = 3.5)
ggarrange(MeHg.vs.DOC, MeHg.vs.SUVA)
dev.off()
