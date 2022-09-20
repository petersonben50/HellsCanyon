#### code/geochem/MeHg_by_redox_conditions.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Read in data ####
geochem.data.adj <- readRDS("dataEdited/geochem/geochem_WC_adjusted_for_redox_classification.rds")


#### Generate needed vectors ####
unique(geochem.data.adj$redox_status)
redox.color.vector <- cb.translator[c("bluishgreen", "reddishpurple", "orange", "skyblue", "blue")]
names(redox.color.vector) <- c("oxic", "suboxic", "no_nitrate_no_sulfide", "no_nitrate_possible_sulfide", "sulfidic")

redox.shape.vector <- c(16, 17, 15, 5, 18)
names(redox.shape.vector) <- names(redox.color.vector)

renaming.vector <- c("Oxygen detected", "No oxygen, nitrate detected", "No nitrate, no sulfide",
                     "No nitrate, sulfide not measured", "Sulfide detected")
names(renaming.vector) <- names(redox.color.vector)


#### Generate plot ####
MeHg.Mn.plot <- geochem.data.adj %>%
  # filter(year(date) == 2019) %>%
  # filter(RM == 286,
  #        elevation_m > 588,
  #        elevation_m < 610) %>%
  # filter(redox_status != "oxic") %>%
  ggplot(aes(x = f_mn_mg_per_l,
             y = MeHg_diss_ngL)) +
  geom_smooth(method = lm ,
              color = "black",
              fill = "grey75",
              se = TRUE,
              level = 0.95) +
  geom_point(aes(color = redox_status,
                 shape = redox_status)) +
  geom_hline(yintercept = 0.01,
             linetype = 2) +
  geom_vline(xintercept = 0.0002,
             linetype = 2) +
  scale_color_manual(values = redox.color.vector,
                     labels = renaming.vector,
                     name = "Redox status") +
  scale_shape_manual(values = redox.shape.vector,
                     labels = renaming.vector,
                     name = "Redox status") +
  theme_classic() +
  scale_x_continuous(limits = c(0.0002, 2),
                     trans = 'log10') +
  scale_y_continuous(limits = c(0.01, 4),
                     trans = 'log10') +
  labs(x = "Filter-passing Mn (mg/L)",
       y = "Dissolved MeHg (ng/L)") +
  theme(legend.position = c(0.38, 0.85),
        legend.box.background = element_rect(colour = "black",
                                             size = 1.2),
        legend.key.size = unit(5, units = "mm"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
MeHg.Mn.plot
# Linearity looks pretty good



#### Linear model with a log-log transformation ####
linear.model.log.log <- lm(log(MeHg_diss_ngL, 10) ~ log(f_mn_mg_per_l, 10),
                           data = geochem.data.adj)
# Check residuals
shapiro.test(linear.model.log.log$residuals)
par(mfrow = c(1,2))
plot(density(linear.model.log.log$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# QQ-normal plot
qqnorm(linear.model.log.log$residuals)
qqline(linear.model.log.log$residuals)
# Pretty close to normal.

# Summarize model
summary(linear.model.log.log)
# Linear fit is great, R^2 = 0.5814



#### Add linear model to plot ####
MeHg.Mn.plot <- MeHg.Mn.plot +
  geom_abline(slope = coef(linear.model.log.log)[[2]],
              intercept = coef(linear.model.log.log)[[1]]) +
  geom_label(x = -0.6, y = 0.5,
             label = paste("Adjusted r2 = ", round(summary(linear.model.log.log)$adj.r.squared, 3), "\n",
                           "p << 0.0001",
                           sep = ""),
             color = "black")
MeHg.Mn.plot


#### Save out plot ####
pdf("results/geochem/MeHgVsMn_redoxStatus.pdf",
    height = 6,
    width = 6)
MeHg.Mn.plot
dev.off()
