#### code/geochem/HgT_time_course.R ####
# Benjamin D. Peterson


#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(ggpubr)
library(lubridate)
library(readxl)
library(tidyverse)
source("code/geochem/profile_functions.R")


#### Read in geochem data ####
geochem.data <- read.csv("dataEdited/geochem/geochem_WC.csv") %>%
  filter(depth != "SW") %>%
  mutate(depth = as.numeric(depth)) %>%
  arrange(depth)
geochem.data.of.interest <- geochem.data %>%
  mutate(RM_date = paste(RM, date, sep = "-")) %>%
  filter(RM_date %in% (geochem.data %>%
                         filter(month(date) >=5 & month(date) <= 11,
                                constituent == "MeHg_diss_percent") %>%
                         group_by(date, RM, constituent) %>%
                         summarise(count = n()) %>%
                         ungroup() %>%
                         filter(count >= 3) %>%
                         mutate(RM_date = paste(RM, date, sep = "-")) %>%
                         select(RM_date) %>%
                         unlist(use.names = FALSE)
                       )) %>%
  select(-RM_date) %>%
  filter(year(date) != 2015)
sampling.dates <- unique(geochem.data.of.interest$date)



#### Generate needed color ramp ####
colorize.function = colorRampPalette(cb.translator[c("black", "orange", "yellow")])
nitrate.vector = colorize.function(7)
names(nitrate.vector) <- floor_date(seq(from = as.Date("2017-05-01"),
                                        to = as.Date("2017-11-01"),
                                        by = 1),
                                    'month') %>%
  unique() %>% month(label = TRUE)


#### Function to generate nitrate plots ####
function.of.function <- function(year.you.want) {
  time.course.profile.plot(geochem.data.to.use = geochem.data.of.interest,
                           years.to.use = year.you.want,
                           concentrations.to.use = c(0, 100),
                           RMs.to.use = c(286,300, 310),
                           parameter.to.plot = "MeHg_diss_percent",
                           xlabel.to.use = "% MeHg",
                           color.vector.to.use = nitrate.vector) +
    facet_wrap(~RM) +
    theme(text = element_text(size = 9))
}
plot.2016 <- function.of.function(2016) + theme(legend.position = "none")
plot.2017 <- function.of.function(2017) + theme(legend.position = "none")
plot.2018 <- function.of.function(2018) + theme(legend.position = "none")
plot.2019 <- function.of.function(2019) + theme(legend.position = "none")


#### Extract legend to plot separately ####
legend.to.use <- get_legend(function.of.function(2016))
# Convert to a ggplot and print
plotted.legend <- as_ggplot(legend.to.use)


#### Save out nitrate plots ####
pdf("results/geochem/MeHg_percent_time_course.pdf",
    width = 7,
    height = 8)
ggarrange(plot.2016, plot.2017, plot.2018, plot.2019,
          plotted.legend,
          labels = c("a", "b", "c", "d"),
          ncol = 2,
          nrow = 3)
dev.off()

