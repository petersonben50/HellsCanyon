#### code/geochem/MeHg_by_DOC_SUVA.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(ggpubr)
library(lubridate)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")
source("code/HCC_plotting_needs.R")


#### Read in data ####
geochem.data.adj <- readRDS("dataEdited/geochem/geochem_WC_adjusted_for_redox_classification.rds") %>%
  filter(year(date) != 2015)


#### Plot DOC by SUVA with redox for overview ####
geochem.data.adj %>%
  filter(!is.na(doc_boulder_mgc_per_l)) %>%
  ggplot(aes(x = doc_boulder_mgc_per_l,
             y = suva_254nm_l_per_mgc_per_m,
             col = redox_status)) +
  geom_point(aes(color = redox_status,
                 shape = as.character(year(date)))) +
  scale_color_manual(values = color.vector,
                     labels = renaming.vector,
                     name = "Redox status") +
  scale_shape_manual(values = shape.vector,
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
  scale_color_manual(values = color.vector,
                     labels = renaming.vector,
                     name = "Redox status") +
  scale_shape_manual(values = shape.vector,
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
  scale_color_manual(values = color.vector,
                     labels = renaming.vector,
                     name = "Redox status") +
  scale_shape_manual(values = shape.vector,
                     labels = renaming.vector,
                     name = "Year") +
  theme_classic() +
  labs(x = "DOC (mg/L)",
       y = "MeHg (ng/L)") +
  theme(legend.position = c(0.3, 0.85))



#### Linear correlation: MeHg to DOC ####
linear.model.doc <- lm(MeHg_diss_ngL ~ doc_boulder_mgc_per_l,
                       data = geochem.data.adj %>%
                         filter(redox_status != "oxic"))
# Summarize model
summary(linear.model.doc)
linear.model.doc <- lm(MeHg_diss_ngL ~ doc_boulder_mgc_per_l,
                       data = geochem.data.adj)
# Summarize model
summary(linear.model.doc)
# Check residuals
# shapiro.test(linear.model.doc$residuals)
# par(mfrow = c(1,2))
plot(density(linear.model.doc$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# QQ-normal plot
qqnorm(linear.model.doc$residuals)
qqline(linear.model.doc$residuals)
# Little bit further from normality, but fairly close.
# Mostly has a right skew here.
# We'll go ahead with it.




#### Linear correlation: MeHg to SUVA ####
linear.model.suva <- lm(MeHg_diss_ngL ~ suva_254nm_l_per_mgc_per_m,
                        data = geochem.data.adj %>%
                          filter(redox_status != "oxic"))
# Summarize model
summary(linear.model.suva)
linear.model.suva <- lm(MeHg_diss_ngL ~ suva_254nm_l_per_mgc_per_m,
                        data = geochem.data.adj)
# Summarize model
summary(linear.model.suva)


# Check residuals
shapiro.test(linear.model.suva$residuals)
par(mfrow = c(1,2))
plot(density(linear.model.suva$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# QQ-normal plot
qqnorm(linear.model.suva$residuals)
qqline(linear.model.suva$residuals)
# Little bit further from normality, but fairly close.
# Mostly has a right skew here.
# We'll go ahead with it.





#### Save out plot ####
pdf("results/geochem/MeHg_DOC.pdf",
    width = 7.2,
    height = 3.5)
ggarrange(MeHg.vs.DOC, MeHg.vs.SUVA)
dev.off()
