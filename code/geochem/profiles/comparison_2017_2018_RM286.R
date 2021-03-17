#### code/profiles/comparison_2017_2018_RM286.R ####
# Benjamin D. Peterson

#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(naniar)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Read in geochem data ####

geochem.data <- read.csv("dataEdited/waterChemistry/geochem_WC_2015_2018.csv",
                         stringsAsFactors = FALSE) %>%
  filter(year(date) %in% c("2017", "2018")) %>%
  filter(month(date) == "9") %>%
  arrange(depth)
# Replace non-detects with zeroes in sulfide column
geochem.data$sulfide_mg.L[geochem.data$sulfide_mg.L == "nd"] <- 0
geochem.data$sulfide_mg.L <- as.numeric(geochem.data$sulfide_mg.L)



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

pdf("results/geochem/profiles/compare_2017_2018_profiles_RM286.pdf",
    height = 6,
    width = 8)

par(mfrow = c(2, 4),
    mgp=c(1.5,0.4,0),
    tck=-0.008,
    mar = c(6, 4.5, 3, 1))
RM.of.interest <- 286
max.depth <- 75

#### First plot: Anions 2017 ####
geochem.data.site.2017 <- geochem.data %>%
  filter(RM == RM.of.interest,
         year(date) == "2017") %>%
  replace_with_na_all(condition = ~.x == "na") %>%
  select(c(depth, Mn.part_ug.L, Mn_ug.L, Fe.part_ug.L, Fe_ug.L, SO4_mg.L, NO3_mgN.L, sulfide_mg.L)) %>%
  as.data.frame()

plot(x = geochem.data.site.2017$NO3_mgN.L,
     y = geochem.data.site.2017$depth,
     xlab = "Nitrate (mgN/L)",
     ylab = "Depth (m)",
     xlim = c(0, 1.2),
     ylim = c(max.depth, 0),
     col = cb.translator["vermillion"],
     pch = 18,
     cex = 1.2,
     main = "Nitrate")
lines(x = geochem.data.site.2017$NO3_mgN.L,
      y = geochem.data.site.2017$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)
title(ylab = paste(2017),
      line = 2.8,
      cex.lab = 2)
par(mar = c(6, 3, 3, 1))

#### Second plot: Fe, 2017 ####
plot(x = geochem.data.site.2017$Fe_ug.L,
     y = geochem.data.site.2017$depth,
     xlab = "Fe (ug/L)",
     ylab = "Depth (m)",
     xlim = c(0, 100),
     ylim = c(max.depth, 0),
     col = cb.translator["orange"],
     pch = 18,
     main = "Fe")
lines(x = geochem.data.site.2017$Fe_ug.L,
      y = geochem.data.site.2017$depth,
      col = cb.translator["orange"],
      lwd = 0.8)
# Add in particulate Fe
points(x = geochem.data.site.2017$Fe.part_ug.L - geochem.data.site.2017$Fe_ug.L,
       y = geochem.data.site.2017$depth,
       col = "orange",
       pch = 5)
lines(x = geochem.data.site.2017$Fe.part_ug.L - geochem.data.site.2017$Fe_ug.L,
      y = geochem.data.site.2017$depth,
      col = cb.translator["orange"],
      lwd = 0.8,
      lty = 2)
# Add legend
legend("bottomright",
       legend = c("Fe (diss.)",
                  "Fe (part.)"),
       col = c("orange",
               "orange"),
       pch = c(18, 5),
       bty = "n")

#### Third plot: Mn, 2017 ####
plot(x = geochem.data.site.2017$Mn_ug.L,
     y = geochem.data.site.2017$depth,
     xlab = "Mn (ug/L)",
     ylab = "Depth (m)",
     xlim = c(0, 1600),
     ylim = c(max.depth, 0),
     pch = 18,
     col = cb.translator["bluishgreen"],
     main = "Mn")
lines(x = geochem.data.site.2017$Mn_ug.L,
      y = geochem.data.site.2017$depth,
      col = cb.translator["bluishgreen"],
      lwd = 0.8)

# Add in particulate Mn
points(x = geochem.data.site.2017$Mn.part_ug.L - geochem.data.site.2017$Mn_ug.L,
       y = geochem.data.site.2017$depth,
       col = cb.translator["bluishgreen"],
       pch = 5)
lines(x = geochem.data.site.2017$Mn.part_ug.L - geochem.data.site.2017$Mn_ug.L,
      y = geochem.data.site.2017$depth,
      col = cb.translator["bluishgreen"],
      lwd = 0.8,
      lty = 2)


# Add legend
legend("topright",
       legend = c("Mn (diss.)",
                  "Mn (part.)"),
       col = c(cb.translator["bluishgreen"],
               cb.translator["bluishgreen"]),
       pch = c(18, 5),
       bty = "n")




#### Fourth plot: Sulfide, 2017 ####
plot(x = geochem.data.site.2017$sulfide_mg.L,
     y = geochem.data.site.2017$depth,
     xlab = "Sulfide (mg/L)",
     ylab = "Depth (m)",
     xlim = c(0, 1.5),
     ylim = c(max.depth, 0),
     pch = 18,
     col = cb.translator["blue"],
     main = "Sulfide")
lines(x = geochem.data.site.2017$sulfide_mg.L,
      y = geochem.data.site.2017$depth,
      col = cb.translator["blue"],
      lwd = 0.8,
      lty = 1)

# Add legend
legend("topright",
       legend = c("Sulfide"),
       col = c("blue",
               "orange"),
       pch = c(18),
       bty = "n")


#### Fifth plot: Nitrate, 2018 ####
par(mar = c(6, 4.5, 3, 1))
geochem.data.site.2018 <- geochem.data %>%
  filter(RM == RM.of.interest,
         year(date) == "2018") %>%
  replace_with_na_all(condition = ~.x == "na") %>%
  select(c(depth, Mn.part_ug.L, Mn_ug.L, Fe.part_ug.L, Fe_ug.L, SO4_mg.L, NO3_mgN.L, sulfide_mg.L)) %>%
  as.data.frame()

plot(x = geochem.data.site.2018$NO3_mgN.L,
     y = geochem.data.site.2018$depth,
     xlab = "Nitrate (mgN/L)",
     ylab = "Depth (m)",
     xlim = c(0, 1.2),
     ylim = c(max.depth, 0),
     col = cb.translator["vermillion"],
     pch = 18,
     cex = 1.2,
     main = "Nitrate")
lines(x = geochem.data.site.2018$NO3_mgN.L,
      y = geochem.data.site.2018$depth,
      col = cb.translator["vermillion"],
      lwd = 0.8)

title(ylab = paste(2018),
      line = 2.8,
      cex.lab = 2)
par(mar = c(6, 3, 3, 1))


#### Sixth plot: Fe, 2018 ####
plot(x = geochem.data.site.2018$Fe_ug.L,
     y = geochem.data.site.2018$depth,
     xlab = "Fe (ug/L)",
     ylab = "Depth (m)",
     xlim = c(0, 100),
     ylim = c(max.depth, 0),
     col = cb.translator["orange"],
     pch = 18,
     main = "Fe")
lines(x = geochem.data.site.2018$Fe_ug.L,
      y = geochem.data.site.2018$depth,
      col = cb.translator["orange"],
      lwd = 0.8)
# Add in particulate Fe
points(x = geochem.data.site.2018$Fe.part_ug.L - geochem.data.site.2018$Fe_ug.L,
       y = geochem.data.site.2018$depth,
       col = "orange",
       pch = 5)
lines(x = geochem.data.site.2018$Fe.part_ug.L - geochem.data.site.2018$Fe_ug.L,
      y = geochem.data.site.2018$depth,
      col = cb.translator["orange"],
      lwd = 0.8,
      lty = 2)
# Add legend
legend("bottomright",
       legend = c("Fe (diss.)",
                  "Fe (part.)"),
       col = c("orange",
               "orange"),
       pch = c(18, 5),
       bty = "n")


#### Seventh plot: Mn, 2018 ####
plot(x = geochem.data.site.2018$Mn_ug.L,
     y = geochem.data.site.2018$depth,
     xlab = "Mn (ug/L)",
     ylab = "Depth (m)",
     xlim = c(0, 1600),
     ylim = c(max.depth, 0),
     pch = 18,
     col = cb.translator["bluishgreen"],
     main = "Mn")
lines(x = geochem.data.site.2018$Mn_ug.L,
      y = geochem.data.site.2018$depth,
      col = cb.translator["bluishgreen"],
      lwd = 0.8)
# Add in particulate Mn
points(x = geochem.data.site.2018$Mn.part_ug.L - geochem.data.site.2018$Mn_ug.L,
       y = geochem.data.site.2018$depth,
       col = cb.translator["bluishgreen"],
       pch = 5)
lines(x = geochem.data.site.2018$Mn.part_ug.L - geochem.data.site.2018$Mn_ug.L,
      y = geochem.data.site.2018$depth,
      col = cb.translator["bluishgreen"],
      lwd = 0.8,
      lty = 2)
# Add legend
legend("topright",
       legend = c("Mn (diss.)",
                  "Mn (part.)"),
       col = c(cb.translator["bluishgreen"],
               cb.translator["bluishgreen"]),
       pch = c(18, 5),
       bty = "n")




#### Eighth plot: sulfide data ####
geochem.data.site.2018.sulfide <- geochem.data.site.2018 %>%
  filter(sulfide_mg.L != "NA")
plot(x = geochem.data.site.2018.sulfide$sulfide_mg.L,
     y = geochem.data.site.2018.sulfide$depth,
     xlab = "Sulfide (mg/L)",
     ylab = "Depth (m)",
     xlim = c(0, 1.5),
     ylim = c(max.depth, 0),
     pch = 18,
     col = cb.translator["blue"],
     main = "Sulfide")
lines(x = geochem.data.site.2018.sulfide$sulfide_mg.L,
      y = geochem.data.site.2018.sulfide$depth,
      col = cb.translator["blue"],
      lwd = 0.8,
      lty = 1)
# Add legend
legend("topright",
       legend = c("Sulfide"),
       col = c("blue",
               "orange"),
       pch = c(18),
       bty = "n")

dev.off()