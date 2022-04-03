#### code/manuscript_figures/supplemental/hg_RM300.R ####
# Benjamin D. Peterson

#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Read in Hg data ####
Hg.data <- read.csv("dataEdited/waterChemistry/Hg_data.csv",
                    stringsAsFactors = FALSE) %>%
  filter(RM == 300,
         year(date) %in% c("2017", "2018"),
         month(date) == "9",
         day(date) > 20) %>%
  group_by(RM, date, depth, constituent) %>%
  summarize(concentration = mean(concentration)) %>%
  arrange(depth) %>%
  spread(key = constituent,
         value = concentration)


#### Read in geochem data ####
geochem.data <- read.csv("dataEdited/waterChemistry/geochem_WC.csv",
                         stringsAsFactors = FALSE) %>%
  filter((year(date) == "2017" & month(date) == "9") |
           (year(date) == "2018" & month(date) == "9" & day(date) > "15")) %>%
  filter(RM == 300) %>%
  arrange(depth) %>%
  group_by(RM, depth, date, constituent) %>%
  summarize(concentration = mean(concentration)) %>%
  spread(key = constituent,
         value = concentration)


#### Read in sonde data ####
seabird.data <- readRDS("dataEdited/seabird/seabird_data.rds") %>%
  filter(year(date) %in% c(2017, 2018),
         RM == 300)



#### Set up PDF with plots ####
pdf("results/manuscript_figures/supplemental_hg/RM300_hg.pdf",
    height = 6,
    width = 5)
par(mfrow = c(2, 4),
    mgp=c(1.5,0.4,0),
    tck=-0.008,
    mar = c(6, 2, 3, .25))
max.depth <- 55


#### 2017 data prep ####
Hg.data.2017 <- Hg.data %>%
  filter(year(date) == "2017") %>%
  replace_with_na_all(condition = ~.x == "na")
geochem.data.2017 <- geochem.data %>%
  filter(year(date) == "2017")
seabird.data.2017 <- seabird.data %>%
  filter(year(date) == "2017")


#### First plot: Redox data, 2017 ####
plot(x = seabird.data.2017$DO,
     y = seabird.data.2017$depth,
     xlab = "DO (mg/L)",
     ylab = "",
     xlim = c(0, 7.5),
     ylim = c(max.depth, 0),
     col = cb.translator["black"],
     pch = 18,
     cex = 1.2,
     main = "")
title(ylab = paste(2017),
      line = 2.8,
      cex.lab = 2)
# Add nitrate data
points(x = geochem.data.2017$NO3_mgN.L*5,
       y = geochem.data.2017$depth,
       col = cb.translator["orange"],
       pch = 18,
       lwd = 0.8)
lines(x = geochem.data.2017$NO3_mgN.L*5,
      y = geochem.data.2017$depth,
      col = cb.translator["orange"],
      lwd = 0.8)
axis(side = 3,
     at = seq(0, 6, by = 1),
     labels = seq(0, 6, by = 1) / 5)

# Add manganese data
points(x = geochem.data.2017$Mn_ug.L/300,
       y = geochem.data.2017$depth,
       col = "purple",
       pch = 18,
       lwd = 0.8)
lines(x = geochem.data.2017$Mn_ug.L/300,
      y = geochem.data.2017$depth,
      col = "purple",
      lwd = 0.8)
# Add in x axes
axis(1,
     line = 3,
     at = seq(0, 6, by = 1),
     labels = seq(0, 6, by = 1) * 0.3)
mtext(1,
      text = "Mn (mg/L)",
      line = 4.5,
      cex = 0.7)
mtext(3,
      text = "Nitrate (mgN/L)",
      line = 1.5,
      cex = 0.7)
legend(x = 3,
       y = 45,
       legend = c("DO",
                  "Nitrate",
                  "Mn"),
       col = c("black", cb.translator["orange"], "purple"),
       pch = 18,
       bty = "n",
       cex = 0.8)




#### Second plot: iHg, 2017 ####
plot(x = Hg.data.2017$FiHg,
     y = Hg.data.2017$depth,
     xlab = "iHg (ng/L)",
     ylab = "",
     xlim = c(0, 1.5),
     ylim = c(max.depth, 0),
     col = cb.translator["bluishgreen"],
     pch = 18,
     cex = 1.2,
     main = "iHg")
lines(x = Hg.data.2017$FiHg,
      y = Hg.data.2017$depth,
      col = cb.translator["bluishgreen"],
      lwd = 0.8)
points(x = Hg.data.2017$PiHg,
       y = Hg.data.2017$depth,
       col = cb.translator["bluishgreen"],
       pch = 5,
       lwd = 0.8)
lines(x = Hg.data.2017$PiHg,
      y = Hg.data.2017$depth,
      col = cb.translator["bluishgreen"],
      lwd = 0.8)
legend("topright",
       legend = c("Dissolved",
                  "Particulate"),
       col = cb.translator["bluishgreen"],
       pch = c(18, 5),
       bty = "n",
       cex = 0.8)


#### Third plot: MeHg, 2017 ####
plot(x = Hg.data.2017$FMHG,
     y = Hg.data.2017$depth,
     xlab = "MeHg (ng/L)",
     ylab = "",
     xlim = c(0, 2.5),
     ylim = c(max.depth, 0),
     col = cb.translator["vermillion"],
     pch = 18,
     cex = 1.2,
     main = "MeHg")
lines(x = Hg.data.2017$FMHG,
      y = Hg.data.2017$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)
points(x = Hg.data.2017$PMHG,
       y = Hg.data.2017$depth,
       col = cb.translator["vermillion"],
       pch = 5,
       lwd = 0.8)
lines(x = Hg.data.2017$PMHG,
      y = Hg.data.2017$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)
legend("topright",
       legend = c("Dissolved",
                  "Particulate"),
       col = cb.translator["vermillion"],
       pch = c(18, 5),
       bty = "n",
       cex = 0.8)


#### Fourth plot: Fraction MeHg, 2017 ####
plot(x = Hg.data.2017$FMHG_per,
     y = Hg.data.2017$depth,
     xlab = "Fraction MeHg",
     ylab = "",
     xlim = c(0, 1),
     ylim = c(max.depth, 0),
     col = cb.translator["black"],
     pch = 18,
     cex = 1.2,
     main = "MeHg")
lines(x = Hg.data.2017$FMHG_per,
      y = Hg.data.2017$depth,
      col = cb.translator["black"],
      lwd = 0.8)
points(x = Hg.data.2017$PMHG_per,
       y = Hg.data.2017$depth,
       col = cb.translator["black"],
       pch = 5,
       lwd = 0.8)
lines(x = Hg.data.2017$PMHG_per,
      y = Hg.data.2017$depth,
      col = cb.translator["black"],
      lwd = 0.8)
# Add legend
legend("topright",
       legend = c("Dissolved",
                  "Particulate"),
       col = "black",
       pch = c(18, 5),
       bty = "n",
       cex = 0.8)





#### 2018 data prep ####
Hg.data.2018 <- Hg.data %>%
  filter(year(date) == "2018") %>%
  replace_with_na_all(condition = ~.x == "na")
geochem.data.2018 <- geochem.data %>%
  filter(year(date) == "2018")
seabird.data.2018 <- seabird.data %>%
  filter(year(date) == "2018")



#### First plot: Redox data, 2018 ####
plot(x = seabird.data.2018$DO,
     y = seabird.data.2018$depth,
     xlab = "DO (mg/L)",
     ylab = "",
     xlim = c(0, 7.5),
     ylim = c(max.depth, 0),
     col = cb.translator["black"],
     pch = 18,
     cex = 1.2,
     main = "")
title(ylab = paste(2018),
      line = 2.8,
      cex.lab = 2)
# Add nitrate data
points(x = geochem.data.2018$NO3_mgN.L*5,
       y = geochem.data.2018$depth,
       col = cb.translator["orange"],
       pch = 18,
       lwd = 0.8)
lines(x = geochem.data.2018$NO3_mgN.L*5,
      y = geochem.data.2018$depth,
      col = cb.translator["orange"],
      lwd = 0.8)
axis(side = 3,
     at = seq(0, 6, by = 1),
     labels = seq(0, 6, by = 1) / 5)

# Add manganese data
points(x = geochem.data.2018$Mn_ug.L/300,
       y = geochem.data.2018$depth,
       col = "purple",
       pch = 18,
       lwd = 0.8)
lines(x = geochem.data.2018$Mn_ug.L/300,
      y = geochem.data.2018$depth,
      col = "purple",
      lwd = 0.8)
# Add in x axes
axis(1,
     line = 3,
     at = seq(0, 6, by = 1),
     labels = seq(0, 6, by = 1) * 0.3)
mtext(1,
      text = "Mn (mg/L)",
      line = 4.5,
      cex = 0.7)
mtext(3,
      text = "Nitrate (mgN/L)",
      line = 1.5,
      cex = 0.7)
legend(x = 3,
       y = 45,
       legend = c("DO",
                  "Nitrate",
                  "Mn"),
       col = c("black", cb.translator["orange"], "purple"),
       pch = 18,
       bty = "n",
       cex = 0.8)


#### Second plot: iHg, 2018 ####
plot(x = Hg.data.2018$FiHg,
     y = Hg.data.2018$depth,
     xlab = "iHg (ng/L)",
     ylab = "",
     xlim = c(0, 1.5),
     ylim = c(max.depth, 0),
     col = cb.translator["bluishgreen"],
     pch = 18,
     cex = 1.2,
     main = "iHg")
lines(x = Hg.data.2018$FiHg,
      y = Hg.data.2018$depth,
      col = cb.translator["bluishgreen"],
      lwd = 0.8)
points(x = Hg.data.2018$PiHg,
       y = Hg.data.2018$depth,
       col = cb.translator["bluishgreen"],
       pch = 5,
       lwd = 0.8)
lines(x = Hg.data.2018$PiHg,
      y = Hg.data.2018$depth,
      col = cb.translator["bluishgreen"],
      lwd = 0.8)
legend("topright",
       legend = c("Dissolved",
                  "Particulate"),
       col = cb.translator["bluishgreen"],
       pch = c(18, 5),
       bty = "n",
       cex = 0.8)


#### Third plot: MeHg, 2018 ####
plot(x = Hg.data.2018$FMHG,
     y = Hg.data.2018$depth,
     xlab = "MeHg (ng/L)",
     ylab = "",
     xlim = c(0, 2.5),
     ylim = c(max.depth, 0),
     col = cb.translator["vermillion"],
     pch = 18,
     cex = 1.2,
     main = "MeHg")
lines(x = Hg.data.2018$FMHG,
      y = Hg.data.2018$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)
points(x = Hg.data.2018$PMHG,
       y = Hg.data.2018$depth,
       col = cb.translator["vermillion"],
       pch = 5,
       lwd = 0.8)
lines(x = Hg.data.2018$PMHG,
      y = Hg.data.2018$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)
legend("topright",
       legend = c("Dissolved",
                  "Particulate"),
       col = cb.translator["vermillion"],
       pch = c(18, 5),
       bty = "n",
       cex = 0.8)


#### Fourth plot: Fraction MeHg, 2018 ####
Hg.data.2018.part <- Hg.data.2018 %>%
  filter(PMHG_per < 1)
plot(x = Hg.data.2018$FMHG_per,
     y = Hg.data.2018$depth,
     xlab = "Fraction MeHg",
     ylab = "",
     xlim = c(0, 1),
     ylim = c(max.depth, 0),
     col = cb.translator["black"],
     pch = 18,
     cex = 1.2,
     main = "MeHg")
lines(x = Hg.data.2018$FMHG_per,
      y = Hg.data.2018$depth,
      col = cb.translator["black"],
      lwd = 0.8)
points(x = Hg.data.2018.part$PMHG_per,
       y = Hg.data.2018.part$depth,
       col = cb.translator["black"],
       pch = 5,
       lwd = 0.8)
lines(x = Hg.data.2018.part$PMHG_per,
      y = Hg.data.2018.part$depth,
      col = cb.translator["black"],
      lwd = 0.8)
# Add legend
legend("topright",
       legend = c("Dissolved",
                  "Particulate"),
       col = "black",
       pch = c(18, 5),
       bty = "n",
       cex = 0.8)




dev.off()