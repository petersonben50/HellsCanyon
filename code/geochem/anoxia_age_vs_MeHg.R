#### code/geochem/anoxia_age_vs_MeHg.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(ggpubr)
library(lubridate)
library(tidyverse)
source("code/HCC_plotting_needs.R")


#### Read in data ####
anoxia.start.df <- readRDS("dataEdited/geochem/anoxia_start_dates.rds") %>%
  filter(sampling.year != 2015)
all.redox.data <- readRDS("dataEdited/geochem/geochem_WC_adjusted_for_redox_classification.rds") %>%
  mutate(sampling.year = as.character(year(date))) %>%
  filter(sampling.year != 2015)


#### Calculate age of anoxia for each sample ####
all.redox.data.year <- full_join(anoxia.start.df,
                                 all.redox.data) %>%
  mutate(anoxia.age = date - date.of.initial.anoxia)


#### Set negative values to zero ####
all.redox.data.year$anoxia.age[all.redox.data.year$anoxia.age < 0] <- 0


#### If a sample is oxygenated, set anoxia age to 0 ####
all.redox.data.year$anoxia.age[all.redox.data.year$diss_oxy_mg_per_l >= 0.5] <- 0


#### Generate plot of anoxia age vs. MeHg ####
plot.MeHg.anoxia <- function(year.of.interest) {
  all.redox.data.year %>%
    filter(sampling.year == year.of.interest) %>%
    ggplot(aes(x = as.integer(anoxia.age),
               y = MeHg_diss_ngL)) +
    geom_point(aes(color = redox_status),
               size = 2) +
    scale_color_manual(values = color.vector,
                       labels = renaming.vector,
                       name = "Redox status") +
    labs(title = year.of.interest,
         x = "Length of anoxia (days)",
         y = "Filter-passing MeHg (ng/L)") +
    theme_classic() +
    scale_y_continuous(limits = c(0.01, 3.5),
                       trans = 'log10') +
    theme(axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12)) 
}


#### Extract legend to plot separately ####
legend.to.use <- get_legend(plot.MeHg.anoxia(2016))
# Convert to a ggplot and print
plotted.legend <- as_ggplot(legend.to.use)


#### Plot graphs ####
plot.2016 <- plot.MeHg.anoxia(2016) + theme(legend.position = "none")
plot.2017 <- plot.MeHg.anoxia(2017) + theme(legend.position = "none")
plot.2018 <- plot.MeHg.anoxia(2018) + theme(legend.position = "none")
plot.2019 <- plot.MeHg.anoxia(2019) + theme(legend.position = "none")


#### Save PDF of graphs ####
pdf("results/geochem/MeHg_age_of_anoxia_log.pdf",
    width = 7.2,
    height = 9)
ggarrange(plot.2016, plot.2017,
          plot.2018, plot.2019,
          plotted.legend,
          nrow = 3,
          ncol = 2,
          labels = c("a.", "b.", "c.", "d."))
dev.off()



#### Statistically test linear correlations between MeHg and anoxia age ####
linear.model <- lm(MeHg_diss_ngL ~ anoxia.age,
                   data = all.redox.data.year %>%
                     filter(redox_status != "oxic"))
summary(linear.model)

linear.model <- lm(MeHg_diss_ngL ~ anoxia.age + sampling.year,
                   data = all.redox.data.year %>%
                     filter(redox_status != "oxic"))
summary(linear.model)
