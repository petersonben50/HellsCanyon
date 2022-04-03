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
  filter(RM %in% c(300, 310),
         year(date) == "2019",
         month(date) == "7",
         day(date) > 20) %>%
  group_by(RM, date, depth, constituent) %>%
  summarize(concentration = mean(concentration)) %>%
  arrange(depth) %>%
  spread(key = constituent,
         value = concentration)


#### Read in geochem data ####
geochem.data <- read.csv("dataEdited/waterChemistry/geochem_WC.csv",
                         stringsAsFactors = FALSE) %>%
  filter(year(date) == "2019" & month(date) == "7") %>%
  filter(RM %in% c(300, 310)) %>%
  arrange(depth) %>%
  group_by(RM, depth, date, constituent) %>%
  summarize(concentration = mean(concentration)) %>%
  spread(key = constituent,
         value = concentration)


#### Read in sonde data ####
seabird.data <- readRDS("dataEdited/seabird/seabird_data.rds") %>%
  filter(year(date) == 2019,
         RM %in% c(300, 310))



#### Set up PDF with plots ####
pdf("results/manuscript_figures/supplemental_hg/2019_hg.pdf",
    height = 6,
    width = 5)
par(mfrow = c(2, 4),
    mgp=c(1.5,0.4,0),
    tck=-0.008,
    mar = c(6, 2, 3, .25))
max.depth <- 58


#### RM300 data prep ####
Hg.data.300 <- Hg.data %>%
  filter(RM == "300") %>%
  replace_with_na_all(condition = ~.x == "na")
geochem.data.300 <- geochem.data %>%
  filter(RM == "300")
seabird.data.300 <- seabird.data %>%
  filter(RM == "300")


#### First plot: Sonde data, 300 ####
plot(x = geochem.data.300$NO3_mgN.L,
     y = geochem.data.300$depth,
     xlab = "Nitrate (mgN/L)",
     ylab = "",
     xlim = c(0, 2),
     ylim = c(max.depth, 0),
     col = cb.translator["orange"],
     pch = 18,
     cex = 1.2,
     main = "")
lines(x = geochem.data.300$NO3_mgN.L,
      y = geochem.data.300$depth,
      col = cb.translator["orange"],
      lwd = 0.8)
# Add manganese data
points(x = geochem.data.300$Mn_ug.L/500,
       y = geochem.data.300$depth,
       col = "purple",
       pch = 18,
       lwd = 0.8)
lines(x = geochem.data.300$Mn_ug.L/500,
      y = geochem.data.300$depth,
      col = "purple",
      lwd = 0.8)
legend(x = 1.1,
       y = 10,
       legend = c("Nitrate",
                  "Mn"),
       col = c(cb.translator["orange"], "purple"),
       pch = 18,
       bty = "n",
       cex = 0.8)
# Add in x axes
axis(1,
     line = 3,
     at = seq(0, 2, by = 0.5),
     labels = seq(0, 2, by = 0.5) * 0.5)
mtext(1,
      text = "Mn (mg/L)",
      line = 4.5,
      cex = 0.7)

#### Second plot: iHg, 300 ####
plot(x = Hg.data.300$FiHg,
     y = Hg.data.300$depth,
     xlab = "iHg (ng/L)",
     ylab = "",
     xlim = c(0, 1),
     ylim = c(max.depth, 0),
     col = cb.translator["bluishgreen"],
     pch = 18,
     cex = 1.2,
     main = "iHg")
lines(x = Hg.data.300$FiHg,
      y = Hg.data.300$depth,
      col = cb.translator["bluishgreen"],
      lwd = 0.8)
points(x = Hg.data.300$PiHg,
       y = Hg.data.300$depth,
       col = cb.translator["bluishgreen"],
       pch = 5,
       lwd = 0.8)
lines(x = Hg.data.300$PiHg,
      y = Hg.data.300$depth,
      col = cb.translator["bluishgreen"],
      lwd = 0.8)
legend("topright",
       legend = c("Dissolved",
                  "Particulate"),
       col = cb.translator["bluishgreen"],
       pch = c(18, 5),
       bty = "n",
       cex = 0.8)


#### Third plot: MeHg, 300 ####
plot(x = Hg.data.300$FMHG,
     y = Hg.data.300$depth,
     xlab = "DO (mg/L)",
     ylab = "",
     xlim = c(0, 1),
     ylim = c(max.depth, 0),
     col = cb.translator["vermillion"],
     pch = 18,
     cex = 1.2,
     main = "MeHg")
lines(x = Hg.data.300$FMHG,
      y = Hg.data.300$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)
points(x = Hg.data.300$PMHG,
       y = Hg.data.300$depth,
       col = cb.translator["vermillion"],
       pch = 5,
       lwd = 0.8)
lines(x = Hg.data.300$PMHG,
      y = Hg.data.300$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)
legend("topright",
       legend = c("Dissolved",
                  "Particulate"),
       col = cb.translator["vermillion"],
       pch = c(18, 5),
       bty = "n",
       cex = 0.8)


#### Fourth plot: Fraction MeHg, 300 ####
plot(x = Hg.data.300$FMHG_per,
     y = Hg.data.300$depth,
     xlab = "Fraction MeHg",
     ylab = "",
     xlim = c(0, 1),
     ylim = c(max.depth, 0),
     col = cb.translator["black"],
     pch = 18,
     cex = 1.2,
     main = "MeHg")
lines(x = Hg.data.300$FMHG_per,
      y = Hg.data.300$depth,
      col = cb.translator["black"],
      lwd = 0.8)
points(x = Hg.data.300$PMHG_per,
       y = Hg.data.300$depth,
       col = cb.translator["black"],
       pch = 5,
       lwd = 0.8)
lines(x = Hg.data.300$PMHG_per,
      y = Hg.data.300$depth,
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





#### RM310 data prep ####
Hg.data.310 <- Hg.data %>%
  filter(RM == "310") %>%
  replace_with_na_all(condition = ~.x == "na")
geochem.data.310 <- geochem.data %>%
  filter(RM == "310")
seabird.data.310 <- seabird.data %>%
  filter(RM == "310")
max.depth <- 40


#### First plot: Sonde data, 310 ####
plot(x = geochem.data.310$NO3_mgN.L,
     y = geochem.data.310$depth,
     xlab = "Nitrate (mgN/L)",
     ylab = "",
     xlim = c(0, 2),
     ylim = c(max.depth, 0),
     col = cb.translator["orange"],
     pch = 18,
     cex = 1.2,
     main = "")
lines(x = geochem.data.310$NO3_mgN.L,
      y = geochem.data.310$depth,
      col = cb.translator["orange"],
      lwd = 0.8)
# Add manganese data
points(x = geochem.data.310$Mn_ug.L/500,
       y = geochem.data.310$depth,
       col = "purple",
       pch = 18,
       lwd = 0.8)
lines(x = geochem.data.310$Mn_ug.L/500,
      y = geochem.data.310$depth,
      col = "purple",
      lwd = 0.8)
legend(x = 1.1,
       y = 10,
       legend = c("Nitrate",
                  "Mn"),
       col = c(cb.translator["orange"], "purple"),
       pch = 18,
       bty = "n",
       cex = 0.8)
# Add in x axes
axis(1,
     line = 3,
     at = seq(0, 2, by = 0.5),
     labels = seq(0, 2, by = 0.5) * 0.5)
mtext(1,
      text = "Mn (mg/L)",
      line = 4.5,
      cex = 0.7)



#### Second plot: iHg, 310 ####
plot(x = Hg.data.310$FiHg,
     y = Hg.data.310$depth,
     xlab = "iHg (ng/L)",
     ylab = "",
     xlim = c(0, 1),
     ylim = c(max.depth, 0),
     col = cb.translator["bluishgreen"],
     pch = 18,
     cex = 1.2,
     main = "iHg")
lines(x = Hg.data.310$FiHg,
      y = Hg.data.310$depth,
      col = cb.translator["bluishgreen"],
      lwd = 0.8)
points(x = Hg.data.310$PiHg,
       y = Hg.data.310$depth,
       col = cb.translator["bluishgreen"],
       pch = 5,
       lwd = 0.8)
lines(x = Hg.data.310$PiHg,
      y = Hg.data.310$depth,
      col = cb.translator["bluishgreen"],
      lwd = 0.8)
legend("topright",
       legend = c("Dissolved",
                  "Particulate"),
       col = cb.translator["bluishgreen"],
       pch = c(18, 5),
       bty = "n",
       cex = 0.8)


#### Third plot: MeHg, 310 ####
plot(x = Hg.data.310$FMHG,
     y = Hg.data.310$depth,
     xlab = "DO (mg/L)",
     ylab = "",
     xlim = c(0, 1),
     ylim = c(max.depth, 0),
     col = cb.translator["vermillion"],
     pch = 18,
     cex = 1.2,
     main = "MeHg")
lines(x = Hg.data.310$FMHG,
      y = Hg.data.310$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)
points(x = Hg.data.310$PMHG,
       y = Hg.data.310$depth,
       col = cb.translator["vermillion"],
       pch = 5,
       lwd = 0.8)
lines(x = Hg.data.310$PMHG,
      y = Hg.data.310$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)
legend("topright",
       legend = c("Dissolved",
                  "Particulate"),
       col = cb.translator["vermillion"],
       pch = c(18, 5),
       bty = "n",
       cex = 0.8)


#### Fourth plot: Fraction MeHg, 310 ####
plot(x = Hg.data.310$FMHG_per,
     y = Hg.data.310$depth,
     xlab = "Fraction MeHg",
     ylab = "",
     xlim = c(0, 1),
     ylim = c(max.depth, 0),
     col = cb.translator["black"],
     pch = 18,
     cex = 1.2,
     main = "MeHg")
lines(x = Hg.data.310$FMHG_per,
      y = Hg.data.310$depth,
      col = cb.translator["black"],
      lwd = 0.8)
points(x = Hg.data.310$PMHG_per,
       y = Hg.data.310$depth,
       col = cb.translator["black"],
       pch = 5,
       lwd = 0.8)
lines(x = Hg.data.310$PMHG_per,
      y = Hg.data.310$depth,
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


