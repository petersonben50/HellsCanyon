#### code/manuscript_figures/supplemental/redox_RM286.R ####
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
  filter((year(date) == "2017" & month(date) == "9") |
           (year(date) == "2018" & month(date) == "9" & day(date) > "15")) %>%
  filter(RM == 286) %>%
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
         year(date) %in% c(2017, 2018),
         RM == 286) %>%
  group_by(RM, depth, date, geneName) %>%
  summarize(coverage = sum(coverage)) %>%
  spread(key = geneName,
         value = coverage)
# EET.data <- readRDS("dataEdited/metabolic_analyses/BBOMP/bbomp_depth_clean.rds")


#### Read in sonde data ####
seabird.data <- readRDS("dataEdited/seabird/seabird_data.rds") %>%
  filter(year(date) %in% c(2017, 2018),
         RM == 286)



# #### Function to generate empty plot ####
# empty.plot <- function(title = "") {
#   plot(x = 0,
#        cex = 0,
#        xaxt = "n",
#        yaxt = "n",
#        xlab = "",
#        ylab = "",
#        bty = "n",
#        main = title)
# }
# 


#### Set up PDF with plots ####
pdf("results/manuscript_figures/supplemental_redox/RM286.pdf",
    height = 5,
    width = 9)
par(mfrow = c(2, 7),
    mgp=c(1.5,0.4,0),
    tck=-0.008,
    mar = c(3, 1.75, 2, .25))
max.depth <- 75


#### 2017 data prep ####
geochem.data.2017 <- geochem.data %>%
  filter(year(date) == "2017") %>%
  replace_with_na_all(condition = ~.x == "na")
teap.data.2017 <- teap.data %>%
  filter(year(date) == "2017")
seabird.data.2017 <- seabird.data %>%
  filter(year(date) == "2017")



#### First plot: Sonde data ####
plot(x = seabird.data.2017$temp,
     y = seabird.data.2017$depth,
     xlab = "Temperature (C)",
     ylab = "",
     xlim = c(0, 25),
     ylim = c(max.depth, 0),
     col = cb.translator["blue"],
     pch = 18,
     cex = 1.2,
     main = "")
points(x = seabird.data.2017$DO*4,
       y = seabird.data.2017$depth,
       col = "black",
       pch = 16)
title(ylab = paste(2017),
      line = 2.8,
      cex.lab = 2)
axis(side = 3,
     at = seq(0, 25, by = 5),
     labels = seq(0, 25, by = 5) / 4)


#### Second plot: N reduction 2017 ####
plot(x = geochem.data.2017$NO3_mgN.L,
     y = geochem.data.2017$depth,
     xlab = "Nitrate (mgN/L)",
     ylab = "",
     xlim = c(0, 1.5),
     ylim = c(max.depth, 0),
     col = cb.translator["vermillion"],
     pch = 18,
     cex = 1.2,
     main = "Nitrate")
lines(x = geochem.data.2017$NO3_mgN.L,
      y = geochem.data.2017$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)

#### Third plot: N reductases 2017 ####
plot(x = teap.data.2017$narG,
     y = teap.data.2017$depth,
     xlab = "Gene abundance",
     ylab = "",
     xlim = c(0, 50),
     ylim = c(max.depth, 0),
     col = cb.translator["vermillion"],
     pch = 18,
     cex = 1.2,
     main = "narG")
lines(x = teap.data.2017$narG,
      y = teap.data.2017$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)


#### Fourth plot: Mn 2017 ####
# Plot dissolved Mn
plot(x = geochem.data.2017$Mn_ug.L,
     y = geochem.data.2017$depth,
     xlab = "Mn (ug/L)",
     ylab = "",
     xlim = c(0, 1600),
     ylim = c(max.depth, 0),
     col = "purple",
     pch = 18,
     cex = 1.2,
     main = "Mn")
lines(x = geochem.data.2017$Mn_ug.L,
      y = geochem.data.2017$depth,
      col = "purple",
      lwd = 0.8)
# Add particulate Mn
points(x = geochem.data.2017$Mn.part_ug.L,
       y = geochem.data.2017$depth,
       col = "purple",
       pch = 5)
lines(x = geochem.data.2017$Mn.part_ug.L,
      y = geochem.data.2017$depth,
      col = "purple",
      lwd = 0.8)
# Add legend
legend("topright",
       legend = c("Mn (diss.)",
                  "Mn (part.)"),
       col = "purple",
       pch = c(18, 5),
       bty = "n")




#### Fifth plot: Fe 2017 ####
# Plot dissolved Fe
plot(x = geochem.data.2017$Fe_ug.L,
     y = geochem.data.2017$depth,
     xlab = "Fe (ug/L)",
     ylab = "",
     xlim = c(0, 100),
     ylim = c(max.depth, 0),
     col = cb.translator["orange"],
     pch = 18,
     cex = 1.2,
     main = "Fe")
lines(x = geochem.data.2017$Fe_ug.L,
      y = geochem.data.2017$depth,
      col = cb.translator["orange"],
      lwd = 0.8)
# Add particulate Fe
points(x = geochem.data.2017$Fe.part_ug.L,
       y = geochem.data.2017$depth,
       col = cb.translator["orange"],
       pch = 5)
lines(x = geochem.data.2017$Fe.part_ug.L,
      y = geochem.data.2017$depth,
      col = cb.translator["orange"],
      lwd = 0.8)
# Add legend
legend(x = 10, y = 50,
       legend = c("Fe (diss.)",
                  "Fe (part.)"),
       col = cb.translator["orange"],
       pch = c(18, 5),
       bty = "n")

# 
# #### Sixth plot: EET genes 2017 ####
# plot(x = geochem.data.2017$NO3_mgN.L,
#      y = geochem.data.2017$depth,
#      cex = 0)
# 

#### Seventh plot: Sulfide ####
plot(x = geochem.data.2017$sulfide_mg.L,
     y = geochem.data.2017$depth,
     xlab = "Sulfide (mg/L)",
     ylab = "",
     xlim = c(0, 1.5),
     ylim = c(max.depth, 0),
     col = cb.translator["blue"],
     pch = 18,
     cex = 1.2,
     main = "Sulfide")
lines(x = geochem.data.2017$sulfide_mg.L,
      y = geochem.data.2017$depth,
      col = cb.translator["blue"],
      lwd = 0.8)


####  Eighth plot: dsrA and dsrD, 2017 ####
plot(x = teap.data.2017$dsrA,
     y = teap.data.2017$depth,
     xlab = "Gene abundance",
     ylab = "",
     xlim = c(0, 4.5),
     ylim = c(max.depth, 0),
     col = cb.translator["blue"],
     pch = 18,
     cex = 1.2,
     main = "Reductive dsr")
lines(x = teap.data.2017$dsrA,
      y = teap.data.2017$depth,
      col = cb.translator["blue"],
      lwd = 0.8)
# Add dsrD
points(x = teap.data.2017$dsrD,
       y = teap.data.2017$depth,
       col = cb.translator["bluishgreen"],
       pch = 17)
lines(x = teap.data.2017$dsrD,
      y = teap.data.2017$depth,
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


#### 2018 data prep ####
max.depth <- 75
geochem.data.2018 <- geochem.data %>%
  filter(year(date) == "2018") %>%
  replace_with_na_all(condition = ~.x == "na")
teap.data.2018 <- teap.data %>%
  filter(year(date) == "2018")
seabird.data.2018 <- seabird.data %>%
  filter(year(date) == "2018")


#### First plot: Sonde data ####
plot(x = seabird.data.2018$temp,
     y = seabird.data.2018$depth,
     xlab = "Temperature ()",
     ylab = "",
     xlim = c(0, 25),
     ylim = c(max.depth, 0),
     col = cb.translator["blue"],
     pch = 18,
     cex = 1.2,
     main = "")
points(x = seabird.data.2018$DO*4,
       y = seabird.data.2018$depth,
       col = "black",
       pch = 16)
axis(side = 3,
     at = seq(0, 25, by = 5),
     labels = seq(0, 25, by = 5) / 4)



#### 10th plot: Nitrate 2018 ####
plot(x = geochem.data.2018$NO3_mgN.L,
     y = geochem.data.2018$depth,
     xlab = "Nitrate (mgN/L)",
     ylab = "",
     xlim = c(0, 1.5),
     ylim = c(max.depth, 0),
     col = cb.translator["vermillion"],
     pch = 18,
     cex = 1.2,
     main = "Nitrate")
lines(x = geochem.data.2018$NO3_mgN.L,
      y = geochem.data.2018$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)


#### 11th plot: N reductases 2018 ####
plot(x = teap.data.2018$narG,
     y = teap.data.2018$depth,
     xlab = "Gene abundance",
     ylab = "",
     xlim = c(0, 50),
     ylim = c(max.depth, 0),
     col = cb.translator["vermillion"],
     pch = 18,
     cex = 1.2,
     main = "narG")
lines(x = teap.data.2018$narG,
      y = teap.data.2018$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)


#### 12th plot: Mn 2018 ####
# Plot dissolved Mn
plot(x = geochem.data.2018$Mn_ug.L,
     y = geochem.data.2018$depth,
     xlab = "Mn (ug/L)",
     ylab = "",
     xlim = c(0, 1600),
     ylim = c(max.depth, 0),
     col = "purple",
     pch = 18,
     cex = 1.2,
     main = "Mn")
lines(x = geochem.data.2018$Mn_ug.L,
      y = geochem.data.2018$depth,
      col = "purple",
      lwd = 0.8)
# Add particulate Mn
points(x = geochem.data.2018$Mn.part_ug.L,
       y = geochem.data.2018$depth,
       col = "purple",
       pch = 5)
lines(x = geochem.data.2018$Mn.part_ug.L,
      y = geochem.data.2018$depth,
      col = "purple",
      lwd = 0.8)
# Add legend
legend("topright",
       legend = c("Mn (diss.)",
                  "Mn (part.)"),
       col = "purple",
       pch = c(18, 5),
       bty = "n")

#### 13th plot: Fe 2018 ####
# Plot dissolved Fe
plot(x = geochem.data.2018$Fe_ug.L,
     y = geochem.data.2018$depth,
     xlab = "Fe (ug/L)",
     ylab = "",
     xlim = c(0, 100),
     ylim = c(max.depth, 0),
     col = cb.translator["orange"],
     pch = 18,
     cex = 1.2,
     main = "Fe")
lines(x = geochem.data.2018$Fe_ug.L,
      y = geochem.data.2018$depth,
      col = cb.translator["orange"],
      lwd = 0.8)
# Add particulate Fe
points(x = geochem.data.2018$Fe.part_ug.L,
       y = geochem.data.2018$depth,
       col = cb.translator["orange"],
       pch = 5)
lines(x = geochem.data.2018$Fe.part_ug.L,
      y = geochem.data.2018$depth,
      col = cb.translator["orange"],
      lwd = 0.8)
# Add legend
legend(x = 10, y = 50,
       legend = c("Fe (diss.)",
                  "Fe (part.)"),
       col = cb.translator["orange"],
       pch = c(18, 5),
       bty = "n")

# #### 14th plot: EET genes 2018 ####
# plot(x = geochem.data.2018$NO3_mgN.L,
#      y = geochem.data.2018$depth,
#      cex = 0)


#### 15th plot: Sulfide ####
plot(x = geochem.data.2018$sulfide_mg.L,
     y = geochem.data.2018$depth,
     xlab = "Sulfide (mg/L)",
     ylab = "",
     xlim = c(0, 1.5),
     ylim = c(max.depth, 0),
     col = cb.translator["blue"],
     pch = 18,
     cex = 1.2,
     main = "Sulfide")
lines(x = geochem.data.2018$sulfide_mg.L,
      y = geochem.data.2018$depth,
      col = cb.translator["blue"],
      lwd = 0.8)


####  16th plot: dsrA and dsrD, 2018 ####
plot(x = teap.data.2018$dsrA,
     y = teap.data.2018$depth,
     xlab = "Gene abundance",
     ylab = "",
     xlim = c(0, 4.5),
     ylim = c(max.depth, 0),
     col = cb.translator["blue"],
     pch = 18,
     cex = 1.2,
     main = "Reductive dsr")
lines(x = teap.data.2018$dsrA,
      y = teap.data.2018$depth,
      col = cb.translator["blue"],
      lwd = 0.8)
# Add dsrD
points(x = teap.data.2018$dsrD,
       y = teap.data.2018$depth,
       col = cb.translator["bluishgreen"],
       pch = 17)
lines(x = teap.data.2018$dsrD,
      y = teap.data.2018$depth,
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
