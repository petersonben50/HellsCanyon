#### code/manuscript_figures/supplemental/redox_RM300.R ####
# Benjamin D. Peterson

#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(naniar)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Read in geochem data ####
geochem.data <- read.csv("dataEdited/waterChemistry/geochem_WC.csv",
                         stringsAsFactors = FALSE) %>%
  filter((year(date) == "2019" & month(date) == "7")) %>%
  filter(RM == 310) %>%
  arrange(depth) %>%
  group_by(RM, depth, date, constituent) %>%
  summarize(concentration = mean(concentration)) %>%
  spread(key = constituent,
         value = concentration)


#### Read in TEAP data ####
genes.of.interest <- c("narG", "narH",
                       "dsrA", "dsrD")
teap.data <- readRDS("dataEdited/metabolic_analyses/depth/metabolicProtein_depth_clean.rds") %>%
  filter(geneName %in% genes.of.interest,
         year(date) %in% c(2019),
         RM == 310) %>%
  group_by(RM, depth, date, geneName) %>%
  summarize(coverage = sum(coverage)) %>%
  spread(key = geneName,
         value = coverage)
# EET.data <- readRDS("dataEdited/metabolic_analyses/BBOMP/bbomp_depth_clean.rds")


#### Read in sonde data ####
seabird.data <- readRDS("dataEdited/seabird/seabird_data.rds") %>%
  filter(year(date) %in% c(2019),
         RM == 310)



#### Function to generate empty plot ####
empty.plot <- function(title = "") {
  plot(x = 0,
       cex = 0,
       xaxt = "n",
       yaxt = "n",
       xlab = "",
       ylab = "",
       bty = "n",
       main = title)
}



#### Set up PDF with plots ####
pdf("results/manuscript_figures/supplemental_redox/RM310.pdf",
    height = 2.5,
    width = 8.75)
par(mfrow = c(1, 7),
    mgp=c(1.5,0.4,0),
    tck=-0.008,
    mar = c(3, 1.75, 2, .25))
max.depth <- 40


#### 2019 data prep ####
geochem.data.2019 <- geochem.data %>%
  filter(year(date) == "2019") %>%
  replace_with_na_all(condition = ~.x == "na")
teap.data.2019 <- teap.data %>%
  filter(year(date) == "2019")
seabird.data.2019 <- seabird.data %>%
  filter(year(date) == "2019")



#### First plot: Sonde data ####
plot(x = seabird.data.2019$temp,
     y = seabird.data.2019$depth,
     xlab = "Temperature (C)",
     ylab = "",
     xlim = c(0, 30),
     ylim = c(max.depth, 0),
     col = cb.translator["blue"],
     pch = 18,
     cex = 1.2,
     main = "")
points(x = seabird.data.2019$DO*4,
       y = seabird.data.2019$depth,
       col = "black",
       pch = 16)
title(ylab = paste(2019),
      line = 2.8,
      cex.lab = 2)
axis(side = 3,
     at = seq(0, 30, by = 5),
     labels = seq(0, 30, by = 5) / 4)


#### Second plot: N reduction 2019 ####
plot(x = geochem.data.2019$NO3_mgN.L,
     y = geochem.data.2019$depth,
     xlab = "Nitrate (mgN/L)",
     ylab = "",
     xlim = c(0, 1.6),
     ylim = c(max.depth, 0),
     col = cb.translator["vermillion"],
     pch = 18,
     cex = 1.2,
     main = "Nitrate")
lines(x = geochem.data.2019$NO3_mgN.L,
      y = geochem.data.2019$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)

#### Third plot: N reductases 2019 ####
plot(x = teap.data.2019$narG,
     y = teap.data.2019$depth,
     xlab = "Gene abundance",
     ylab = "",
     xlim = c(0, 60),
     ylim = c(max.depth, 0),
     col = cb.translator["vermillion"],
     pch = 18,
     cex = 1.2,
     main = "narG")
lines(x = teap.data.2019$narG,
      y = teap.data.2019$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)


#### Fourth plot: Mn 2019 ####
# Plot dissolved Mn
plot(x = geochem.data.2019$Mn_ug.L,
     y = geochem.data.2019$depth,
     xlab = "Mn (ug/L)",
     ylab = "",
     xlim = c(0, 1600),
     ylim = c(max.depth, 0),
     col = "purple",
     pch = 18,
     cex = 1.2,
     main = "Mn")
lines(x = geochem.data.2019$Mn_ug.L,
      y = geochem.data.2019$depth,
      col = "purple",
      lwd = 0.8)
# Add particulate Mn
points(x = geochem.data.2019$Mn.part_ug.L,
       y = geochem.data.2019$depth,
       col = "purple",
       pch = 5)
lines(x = geochem.data.2019$Mn.part_ug.L,
      y = geochem.data.2019$depth,
      col = "purple",
      lwd = 0.8)
# Add legend
legend("topright",
       legend = c("Mn (diss.)",
                  "Mn (part.)"),
       col = "purple",
       pch = c(18, 5),
       bty = "n")




#### Fifth plot: Fe 2019 ####
# Plot dissolved Fe
plot(x = geochem.data.2019$Fe_ug.L,
     y = geochem.data.2019$depth,
     xlab = "Fe (ug/L)",
     ylab = "",
     xlim = c(0, 280),
     ylim = c(max.depth, 0),
     col = cb.translator["orange"],
     pch = 18,
     cex = 1.2,
     main = "Fe")
lines(x = geochem.data.2019$Fe_ug.L,
      y = geochem.data.2019$depth,
      col = cb.translator["orange"],
      lwd = 0.8)
# Add particulate Fe
points(x = geochem.data.2019$Fe.part_ug.L,
       y = geochem.data.2019$depth,
       col = cb.translator["orange"],
       pch = 5)
lines(x = geochem.data.2019$Fe.part_ug.L,
      y = geochem.data.2019$depth,
      col = cb.translator["orange"],
      lwd = 0.8)
# Add legend
legend(x = 30, y = 0,
       legend = c("Fe (diss.)",
                  "Fe (part.)"),
       col = cb.translator["orange"],
       pch = c(18, 5),
       bty = "n")


#### Sixth plot: EET genes 2019 ####
plot(x = geochem.data.2019$NO3_mgN.L,
     y = geochem.data.2019$depth,
     cex = 0)


####  Seventh plot: dsrA and dsrD, 2019 ####
plot(x = teap.data.2019$dsrA,
     y = teap.data.2019$depth,
     xlab = "Gene abundance",
     ylab = "",
     xlim = c(0, 4.5),
     ylim = c(max.depth, 0),
     col = cb.translator["blue"],
     pch = 18,
     cex = 1.2,
     main = "Reductive dsr")
lines(x = teap.data.2019$dsrA,
      y = teap.data.2019$depth,
      col = cb.translator["blue"],
      lwd = 0.8)
# Add dsrD
points(x = teap.data.2019$dsrD,
       y = teap.data.2019$depth,
       col = cb.translator["bluishgreen"],
       pch = 17)
lines(x = teap.data.2019$dsrD,
      y = teap.data.2019$depth,
      col = cb.translator["bluishgreen"],
      lwd = 0.8)
# Add legend
legend("topright",
       legend = c("dsrA",
                  "dsrD"),
       col = c(cb.translator["blue"],
               cb.translator["bluishgreen"]),
       pch = c(18, 17),
       bty = "n")


dev.off()
