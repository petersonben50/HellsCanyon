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
         year(date) %in% c(2019)) %>%
  group_by(RM, depth, date, geneName) %>%
  summarize(coverage = sum(coverage)) %>%
  spread(key = geneName,
         value = coverage)

EET.data <- readRDS("dataEdited/metabolic_analyses/BBOMP/bbomp_depth_clean.rds") %>%
  filter(year(date) == 2019,
         classification == "cluster_1_ExtE") %>%
  group_by(RM, depth, date) %>%
  summarise(coverage = sum(coverage)) %>%
  ungroup()


#### Read in sonde data ####
seabird.data <- readRDS("dataEdited/seabird/seabird_data.rds") %>%
  filter(year(date) %in% c(2019))



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
pdf("results/manuscript_figures/supplemental_redox/2019_redox.pdf",
    height = 5,
    width = 8.75)
par(mfrow = c(2, 7),
    mgp=c(1.5,0.4,0),
    tck=-0.008,
    mar = c(3, 1.75, 2, .25))
max.depth <- 60


#### 2019 data prep ####
geochem.data.300 <- geochem.data %>%
  filter(RM == 300) %>%
  replace_with_na_all(condition = ~.x == "na")
teap.data.300 <- teap.data %>%
  filter(RM == 300)
seabird.data.300 <- seabird.data %>%
  filter(RM == 300)
EET.data.300 <- EET.data %>%
  filter(RM == 300)


#### First plot: Sonde data ####
empty.plot()
# plot(x = seabird.data.300$temp,
#      y = seabird.data.300$depth,
#      xlab = "Temperature (C)",
#      ylab = "",
#      xlim = c(0, 30),
#      ylim = c(max.depth, 0),
#      col = cb.translator["blue"],
#      pch = 18,
#      cex = 1.2,
#      main = "")
# points(x = seabird.data.300$DO*4,
#        y = seabird.data.300$depth,
#        col = "black",
#        pch = 16)
# title(ylab = paste(2019),
#       line = 2.8,
#       cex.lab = 2)
# axis(side = 3,
#      at = seq(0, 30, by = 5),
#      labels = seq(0, 30, by = 5) / 4)


#### Second plot: N reduction 2019 ####
plot(x = geochem.data.300$NO3_mgN.L,
     y = geochem.data.300$depth,
     xlab = "Nitrate (mgN/L)",
     ylab = "",
     xlim = c(0, 1.6),
     ylim = c(max.depth, 0),
     col = cb.translator["vermillion"],
     pch = 18,
     cex = 1.2,
     main = "Nitrate")
lines(x = geochem.data.300$NO3_mgN.L,
      y = geochem.data.300$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)

#### Third plot: N reductases 2019 ####
plot(x = teap.data.300$narG,
     y = teap.data.300$depth,
     xlab = "Gene abundance",
     ylab = "",
     xlim = c(0, 60),
     ylim = c(max.depth, 0),
     col = cb.translator["vermillion"],
     pch = 18,
     cex = 1.2,
     main = "narG")
lines(x = teap.data.300$narG,
      y = teap.data.300$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)


#### Fourth plot: Mn 2019 ####
# Plot dissolved Mn
plot(x = geochem.data.300$Mn_ug.L,
     y = geochem.data.300$depth,
     xlab = "Mn (ug/L)",
     ylab = "",
     xlim = c(0, 1600),
     ylim = c(max.depth, 0),
     col = "purple",
     pch = 18,
     cex = 1.2,
     main = "Mn")
lines(x = geochem.data.300$Mn_ug.L,
      y = geochem.data.300$depth,
      col = "purple",
      lwd = 0.8)
# Add particulate Mn
points(x = geochem.data.300$Mn.part_ug.L,
       y = geochem.data.300$depth,
       col = "purple",
       pch = 5)
lines(x = geochem.data.300$Mn.part_ug.L,
      y = geochem.data.300$depth,
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
plot(x = geochem.data.300$Fe_ug.L,
     y = geochem.data.300$depth,
     xlab = "Fe (ug/L)",
     ylab = "",
     xlim = c(0, 280),
     ylim = c(max.depth, 0),
     col = cb.translator["orange"],
     pch = 18,
     cex = 1.2,
     main = "Fe")
lines(x = geochem.data.300$Fe_ug.L,
      y = geochem.data.300$depth,
      col = cb.translator["orange"],
      lwd = 0.8)
# Add particulate Fe
points(x = geochem.data.300$Fe.part_ug.L,
       y = geochem.data.300$depth,
       col = cb.translator["orange"],
       pch = 5)
lines(x = geochem.data.300$Fe.part_ug.L,
      y = geochem.data.300$depth,
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

plot(x = EET.data.300$coverage,
     y = EET.data.300$depth,
     xlab = "Gene abundance",
     ylab = "",
     xlim = c(0, 2),
     ylim = c(max.depth, 0),
     col = cb.translator["orange"],
     pch = 18,
     cex = 1.2,
     main = "Reductive dsr")
lines(x = EET.data.300$coverage,
      y = EET.data.300$depth,
      col = cb.translator["orange"],
      lwd = 0.8)


####  Seventh plot: dsrA and dsrD, 2019 ####
plot(x = teap.data.300$dsrA,
     y = teap.data.300$depth,
     xlab = "Gene abundance",
     ylab = "",
     xlim = c(0, 4.5),
     ylim = c(max.depth, 0),
     col = cb.translator["blue"],
     pch = 18,
     cex = 1.2,
     main = "Reductive dsr")
lines(x = teap.data.300$dsrA,
      y = teap.data.300$depth,
      col = cb.translator["blue"],
      lwd = 0.8)
# Add dsrD
points(x = teap.data.300$dsrD,
       y = teap.data.300$depth,
       col = cb.translator["bluishgreen"],
       pch = 17)
lines(x = teap.data.300$dsrD,
      y = teap.data.300$depth,
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








#### 2019 data prep ####
geochem.data.310 <- geochem.data %>%
  filter(RM == 310) %>%
  replace_with_na_all(condition = ~.x == "na")
teap.data.310 <- teap.data %>%
  filter(RM == 310)
seabird.data.310 <- seabird.data %>%
  filter(RM == 310)
EET.data.310 <- EET.data %>%
  filter(RM == 310)
max.depth <- 40


#### First plot: Sonde data ####
plot(x = seabird.data.310$temp,
     y = seabird.data.310$depth,
     xlab = "Temperature (C)",
     ylab = "",
     xlim = c(0, 30),
     ylim = c(max.depth, 0),
     col = cb.translator["blue"],
     pch = 18,
     cex = 1.2,
     main = "")
points(x = seabird.data.310$DO*4,
       y = seabird.data.310$depth,
       col = "black",
       pch = 16)
title(ylab = paste(2019),
      line = 2.8,
      cex.lab = 2)
axis(side = 3,
     at = seq(0, 30, by = 5),
     labels = seq(0, 30, by = 5) / 4)


#### Second plot: N reduction 2019 ####
plot(x = geochem.data.310$NO3_mgN.L,
     y = geochem.data.310$depth,
     xlab = "Nitrate (mgN/L)",
     ylab = "",
     xlim = c(0, 1.6),
     ylim = c(max.depth, 0),
     col = cb.translator["vermillion"],
     pch = 18,
     cex = 1.2,
     main = "Nitrate")
lines(x = geochem.data.310$NO3_mgN.L,
      y = geochem.data.310$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)

#### Third plot: N reductases 2019 ####
plot(x = teap.data.310$narG,
     y = teap.data.310$depth,
     xlab = "Gene abundance",
     ylab = "",
     xlim = c(0, 60),
     ylim = c(max.depth, 0),
     col = cb.translator["vermillion"],
     pch = 18,
     cex = 1.2,
     main = "narG")
lines(x = teap.data.310$narG,
      y = teap.data.310$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)


#### Fourth plot: Mn 2019 ####
# Plot dissolved Mn
plot(x = geochem.data.310$Mn_ug.L,
     y = geochem.data.310$depth,
     xlab = "Mn (ug/L)",
     ylab = "",
     xlim = c(0, 1600),
     ylim = c(max.depth, 0),
     col = "purple",
     pch = 18,
     cex = 1.2,
     main = "Mn")
lines(x = geochem.data.310$Mn_ug.L,
      y = geochem.data.310$depth,
      col = "purple",
      lwd = 0.8)
# Add particulate Mn
points(x = geochem.data.310$Mn.part_ug.L,
       y = geochem.data.310$depth,
       col = "purple",
       pch = 5)
lines(x = geochem.data.310$Mn.part_ug.L,
      y = geochem.data.310$depth,
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
plot(x = geochem.data.310$Fe_ug.L,
     y = geochem.data.310$depth,
     xlab = "Fe (ug/L)",
     ylab = "",
     xlim = c(0, 280),
     ylim = c(max.depth, 0),
     col = cb.translator["orange"],
     pch = 18,
     cex = 1.2,
     main = "Fe")
lines(x = geochem.data.310$Fe_ug.L,
      y = geochem.data.310$depth,
      col = cb.translator["orange"],
      lwd = 0.8)
# Add particulate Fe
points(x = geochem.data.310$Fe.part_ug.L,
       y = geochem.data.310$depth,
       col = cb.translator["orange"],
       pch = 5)
lines(x = geochem.data.310$Fe.part_ug.L,
      y = geochem.data.310$depth,
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

plot(x = EET.data.310$coverage,
     y = EET.data.310$depth,
     xlab = "Gene abundance",
     ylab = "",
     xlim = c(0, 2),
     ylim = c(max.depth, 0),
     col = cb.translator["orange"],
     pch = 18,
     cex = 1.2,
     main = "Reductive dsr")
lines(x = EET.data.310$coverage,
      y = EET.data.310$depth,
      col = cb.translator["orange"],
      lwd = 0.8)


####  Seventh plot: dsrA and dsrD, 2019 ####
plot(x = teap.data.310$dsrA,
     y = teap.data.310$depth,
     xlab = "Gene abundance",
     ylab = "",
     xlim = c(0, 4.5),
     ylim = c(max.depth, 0),
     col = cb.translator["blue"],
     pch = 18,
     cex = 1.2,
     main = "Reductive dsr")
lines(x = teap.data.310$dsrA,
      y = teap.data.310$depth,
      col = cb.translator["blue"],
      lwd = 0.8)
# Add dsrD
points(x = teap.data.310$dsrD,
       y = teap.data.310$depth,
       col = cb.translator["bluishgreen"],
       pch = 17)
lines(x = teap.data.310$dsrD,
      y = teap.data.310$depth,
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
