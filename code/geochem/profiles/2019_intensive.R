#### code/profiles/2019_intensive.R ####
# Benjamin D. Peterson


#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(naniar)
library(readxl)
library(tidyverse)



#### Read in and prep SeaCat data ####

sonde.data <- read_xlsx("dataRaw/SeaCat/2019_intensive_seacat_brownlee.xlsx") %>%
  select(-c(Datetime, Date, Time,
            `p H`, Location))
names(sonde.data) <- c("depth", "temp", "redox_pot", "cDOM",
                       "Chla", "turb", "PAR", "spec_cond", "DO", "RM")
sonde.data <- sonde.data %>%
  filter(RM != "322")



#### Read in geochem data ####

geochem.data <- read.csv("dataEdited/waterChemistry/geochem_WC.csv",
                         stringsAsFactors = FALSE) %>%
  filter(year(date) == 2019)
geochem.data.OL <- read.csv("dataEdited/waterChemistry/geochem_OL_2019.csv",
                            stringsAsFactors = FALSE)



#### Read in Hg data ####

Hg.data <- read.csv("dataEdited/waterChemistry/Hg_data.csv",
                    stringsAsFactors = FALSE) %>%
  filter(year(date) == 2019)
Hg.data.OL <- read.csv("dataEdited/waterChemistry/Hg_OL_2019_intensive.csv",
                       stringsAsFactors = FALSE)


#### Read in NA sampling information ####
MG.metadata <- read.csv("metadata/metagenome_metadata.csv",
                        stringsAsFactors = FALSE)


#### Read in hgcA depth, add to NA info ####
hgcA.depth <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv",
                       stringsAsFactors = FALSE) %>%
  filter(rep) %>%
  group_by(metagenomeID) %>%
  summarise(coverage = sum(coverage))
hgcA.depth.vector <- hgcA.depth$coverage
names(hgcA.depth.vector) <- hgcA.depth$metagenomeID

MG.depth <- MG.metadata %>%
  mutate(hgcA_depth = hgcA.depth.vector[metagenomeID]) %>%
  filter(year(date) == 2019)


#### List of sites ####
site.list <- sort(unique(c(sonde.data$RM,
                           geochem.data$RM,
                           Hg.data$RM)))


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


#### Plot temp/DO, redox, and MeHg in separate plots data for all sites ####

pdf("results/geochem/profiles/2019_summary_profiles.pdf",
    height = 20,
    width = 15)

par(mfrow = c(6, 8),
    mgp=c(1.5,0.4,0),
    tck=-0.008)

for (RM.of.interest in site.list) {
  
  par(mar = c(6, 4.5, 3, 1))
  

  #### First plot: Temp and DO ####
  if (RM.of.interest %in% sonde.data$RM) {
    
    sonde.data.site <- sonde.data %>%
      filter(RM == RM.of.interest) %>%
      select(-c(RM)) %>%
      as.data.frame()
    
    max.depth <- max(sonde.data.site$depth)
    
    plot(x = sonde.data.site$temp/2.5,
         y = sonde.data.site$depth,
         xaxt = "n",
         xlim = c(0, 10),
         ylim = c(max.depth, 0),
         ylab = "Depth (m)",
         xlab = "Temp (C)",
         pch = 16,
         col = "blue",
         cex = 0.8,
         main = paste("Sonde cast",
                      sep = ""))
    points(x = sonde.data.site$DO,
           y = sonde.data.site$depth,
           pch = 18)
    
    # Add axes
    axis(1,
         at = seq(0, 10, by = 2),
         labels = seq(0, 10, by = 2)*2.5)
    axis(1,
         line = 3,
         at = seq(0, 10, by = 2),
         labels = seq(0, 10, by = 2))
    
    # Add label for DO
    title(xlab = "Dissolved Oxygen (mg/L)",
          line = 4.5)
    
    # Add DNA sampling information
    
    MGs.site <- MG.depth %>%
      filter(RM == RM.of.interest)
    
    if (length(MGs.site$depth) > 0) {
      
      points(x = rep(10, length(MGs.site$depth)),
             y = MGs.site$depth,
             pch = 16,
             col = "darkgreen")
      text(x = rep(6, length(MGs.site$depth)),
           y = MGs.site$depth,
           labels = MGs.site$metagenomeID,
           col = "darkgreen")
    }
    
    
    
    
    # Add legend
    legend("topleft",
           legend = c("Temp (C)",
                      "DO (mg/L)",
                      "NA sample",
                      "MG sequenced"),
           col = c("blue",
                   "black",
                   "darkgreen",
                   "darkgreen"),
           pch = c(16, 16, 1, 16),
           bty = "n")
    
    # Clean up
    rm(sonde.data.site)
    
  } else {
    
    # Get max depth from geochem
    max.depth <- geochem.data %>%
      filter(RM == RM.of.interest) %>%
      select(depth) %>%
      unlist(use.names = FALSE) %>%
      max() %>%
      as.numeric()
    # Plot out empty graph
    plot(x = NA,
         y = NA,
         xaxt = "n",
         xlim = c(0, 10),
         ylim = c(max.depth, 0),
         ylab = "Depth (m)",
         xlab = "",
         pch = 16,
         cex = 0.8,
         main = "Missing sonde cast")
    
    # Plot out the MG data
    
    MGs.site <- MG.depth %>%
      filter(RM == RM.of.interest)
    
    if (length(MGs.site$depth) > 0) {
      
      points(x = rep(10, length(MGs.site$depth)),
             y = MGs.site$depth,
             pch = 16,
             col = "darkgreen")
      text(x = rep(6, length(MGs.site$depth)),
           y = MGs.site$depth,
           labels = MGs.site$metagenomeID,
           col = "darkgreen")
    } else {
      empty.plot()
    }
    
  }
  
  title(ylab = paste("RM",
                     RM.of.interest,
                     sep = ""),
        line = 3,
        cex.lab = 2)
  
  par(mar = c(6, 3, 3, 1))
  
  #### Second plot: Anions ####
  
  if (RM.of.interest %in% geochem.data$RM) {
    
    geochem.data.site <- geochem.data %>%
      filter(RM == RM.of.interest) %>%
      replace_with_na_all(condition = ~.x == "na") %>%
      group_by(RM, depth, date, constituent) %>%
      summarise(concentration = mean(concentration)) %>%
      spread(key = constituent,
             value = concentration) %>%
      as.data.frame() %>%
      ungroup()
    
    plot(x = geochem.data.site$SO4_mg.L/5,
         y = geochem.data.site$depth,
         xlab = "Sulfate (mg/L), Nitrite (ug/L)",
         ylab = "Depth (m)",
         xlim = c(0, 10),
         ylim = c(max.depth, 0),
         xaxt = "n",
         col = "blue",
         pch = 16,
         cex = 1.2,
         main = "Anion concentrations")
    points(x = geochem.data.site$NO3_mgN.L*4,
           y = geochem.data.site$depth,
           col = "red",
           pch = 16)
    points(x = geochem.data.site$NO2_mgN.L*1000/5,
           y = geochem.data.site$depth,
           col = "red",
           pch = 1)
    
    # Add axes
    axis(1,
         at = seq(0, 10, by = 2),
         labels = seq(0, 10, by = 2)*5)
    axis(1,
         line = 3,
         at = seq(0, 10, by = 2),
         labels = seq(0, 10, by = 2)/4)
    
    # Add label for nitrate
    title(xlab = "Nitrate (mgN/L)",
          line = 4.5)
    
    
    # Add legend
    legend("topleft",
           legend = c("Sulfate",
                      "Nitrate",
                      "Nitrite"),
           col = c("blue",
                   "red",
                   "red"),
           pch = c(16, 16, 1),
           cex = 0.9,
           bty = "n")
    
    # Clean up

  } else {
    empty.plot("No anion data")
  }
  
  
  #### Third plot: Anion ratios ####
  if (RM.of.interest %in% geochem.data$RM) {
    
    geochem.data.site <- geochem.data.site %>%
      mutate(SO4.Cl.ratio = SO4_mg.L / Cl_mg.L) %>%
      mutate(NO3.Cl.ratio = NO3_mgN.L / Cl_mg.L) %>%
      as.data.frame()
    
    plot(x = geochem.data.site$SO4.Cl.ratio,
         y = geochem.data.site$depth,
         xlab = "Sulfate/Chloride",
         ylab = "Depth (m)",
         xlim = c(0, 10),
         ylim = c(max.depth, 0),
         col = "blue",
         pch = 16,
         cex = 1.2,
         main = "Anion ratios")
    points(x = geochem.data.site$NO3.Cl.ratio*100,
           y = geochem.data.site$depth,
           col = "red",
           pch = 16)
    
    # Add axes for NO3/Cl ratio
    axis(1,
         line = 3,
         at = seq(0, 10, by = 2),
         labels = seq(0, 10, by = 2)/100)
    
    # Add label for nitrate/chloride ratio
    title(xlab = "Nitrate/Chloride",
          line = 4.5)
    
    # Add legend
    legend("topleft",
           legend = c("Sulfate/Cl-",
                      "Nitrate/Cl-"),
           col = c("blue",
                   "red"),
           pch = c(16, 16, 1),
           cex = 0.9,
           bty = "n")
    
  } else {
    empty.plot("No anion data")
  }
  
  
  #### Fourth and fifth plots: Metal data ####
  if (RM.of.interest %in% geochem.data$RM) {
    # 
    # geochem.data.site <- geochem.data %>%
    #   filter(RM == RM.of.interest) %>%
    #   replace_with_na_all(condition = ~.x == "na") %>%
    #   select(c(depth, Fe_ug.L, Mn_ug.L, Fe.part_ug.L,
    #            Mn.part_ug.L)) %>%
    #   as.data.frame()
    
    geochem.data.OL.site <- geochem.data.OL %>%
      filter(RM == RM.of.interest) %>%
      replace_with_na_all(condition = ~.x == "na") %>%
      select(c(depth, Fe_ug.L, Mn_ug.L, Fe.part_ug.L,
               Mn.part_ug.L)) %>%
      arrange(depth) %>%
      as.data.frame()
    
    
    
    iron.max <- max(c(geochem.data.OL.site$Fe.part_ug.L,
                      geochem.data.OL.site$Fe_ug.L,
                      geochem.data.site$Fe_ug.L,
                      geochem.data.site$Fe.part_ug.L),
                    na.rm = TRUE)
    
    plot(x = geochem.data.site$Fe_ug.L,
         y = geochem.data.site$depth,
         xlab = "Fe (ug/L)",
         ylab = "Depth (m)",
         xlim = c(0, iron.max*1.1),
         ylim = c(max.depth, 0),
         pch = 16,
         cex = 0,
         main = "Iron concentrations")
    
    # Add in dissolved Fe
    points(x = geochem.data.site$Fe_ug.L,
           y = geochem.data.site$depth,
           col = "orange",
           pch = 16)
    points(x = geochem.data.OL.site$Fe_ug.L,
           y = geochem.data.OL.site$depth,
           col = "orange",
           pch = c(16, 17))
    # Add in particulate Fe
    points(x = geochem.data.site$Fe.part_ug.L - geochem.data.site$Fe_ug.L,
           y = geochem.data.site$depth,
           col = "orange",
           pch = 1)
    points(x = geochem.data.OL.site$Fe.part_ug.L - geochem.data.OL.site$Fe_ug.L,
           y = geochem.data.OL.site$depth,
           col = "orange",
           pch = c(1, 2))
    legend("topright",
           legend = c("Fe (diss.)",
                      "Fe (part.)",
                      "Lowest point"),
           col = c("orange",
                   "orange",
                   "black"),
           pch = c(16, 1, 17),
           bty = "n")
    
    
    
    
    Mn.max <- max(c(geochem.data.OL.site$Mn.part_ug.L,
                    geochem.data.OL.site$Mn_ug.L,
                    geochem.data.site$Mn_ug.L,
                    geochem.data.site$Mn.part_ug.L),
                  na.rm = TRUE)
    
    plot(x = geochem.data.site$Mn_ug.L,
         y = geochem.data.site$depth,
         xlab = "Manganese (ug/L)",
         ylab = "Depth (m)",
         xlim = c(0, Mn.max*1.1),
         ylim = c(max.depth, 0),
         pch = 16,
         cex = 0,
         main = "Mn concentrations")
    
    # Add in particulate Mn
    points(x = geochem.data.site$Mn.part_ug.L - geochem.data.site$Mn_ug.L,
           y = geochem.data.site$depth,
           col = "purple",
           pch = 1)
    points(x = geochem.data.OL.site$Mn.part_ug.L - geochem.data.OL.site$Mn_ug.L,
           y = geochem.data.OL.site$depth,
           col = "purple",
           pch = c(1,2))
    # Add in dissolved Mn
    points(x = geochem.data.site$Mn_ug.L,
           y = geochem.data.site$depth,
           col = "purple",
           pch = 16)
    points(x = geochem.data.OL.site$Mn_ug.L,
           y = geochem.data.OL.site$depth,
           col = "purple",
           pch = c(16, 17))
    
    # Add legend
    legend("topright",
           legend = c("Mn (diss.)",
                      "Mn (part.)",
                      "Lowest point"),
           col = c("purple",
                   "purple",
                   "black"),
           pch = c(16, 1, 17),
           bty = "n")
    
    # Clean up
    rm(geochem.data.site)   
    
  } else {
    empty.plot("No anion data")
  }
  
  
  #### 6th and 7th plots: MeHg and iHg ####
  if (RM.of.interest %in% Hg.data$RM) {
    
    Hg.data.site <- Hg.data %>%
      filter(RM == RM.of.interest)
    Hg.data.OL.site <- Hg.data.OL %>%
      filter(RM == RM.of.interest) %>%
      arrange(depth)
    
    plot(x = Hg.data.site$FMHG,
         y = Hg.data.site$depth,
         xlab = "Hg (ng/L)", 
         ylab = "Depth (m)",
         xlim = c(0, 0.5),
         ylim = c(max.depth, 0),
         main = paste("MeHg data", 
                      sep = ""),
         col = "firebrick3",
         pch = 16)
    points(x = Hg.data.OL.site$FMHG,
           y = Hg.data.OL.site$depth,
           col = "firebrick3",
           pch = c(16, 17))
    
    points(x = Hg.data.site$PMHG,
           y = Hg.data.site$depth,
           col = "firebrick3",
           pch = 1)
    points(x = Hg.data.OL.site$PMHG,
           y = Hg.data.OL.site$depth,
           col = "firebrick3",
           pch = c(1, 2))
    
    legend("topright",
           legend = c("Diss. MeHg",
                      "Part. MeHg",
                      "Lowest point"),
           col = c("firebrick3",
                   "firebrick3",
                   "black"),
           pch = c(16,
                   1,
                   17),
           cex = 0.9,
           bty = "n")
    
    # Inorganic Hg
    plot(x = Hg.data.site$FTHG - Hg.data.site$FMHG,
         y = Hg.data.site$depth,
         xlab = "iHg (ng/L)", 
         ylab = "Depth (m)",
         xlim = c(0, 2),
         ylim = c(max.depth, 0),
         main = paste("iHg data", 
                      sep = ""),
         col = "black",
         pch = 16)
    
    points(x = Hg.data.OL.site$FTHG - Hg.data.OL.site$FMHG,
           y = Hg.data.OL.site$depth,
           col = "black",
           pch = 16)
    
    points(x = Hg.data.site$PTHG - Hg.data.site$PMHG,
           y = Hg.data.site$depth,
           col = "black",
           pch = 1)
    points(x = Hg.data.OL.site$PTHG - Hg.data.OL.site$PMHG,
           y = Hg.data.OL.site$depth,
           col = "black",
           pch = 1)
    
    
    legend("topright",
           legend = c("Diss. iHg",
                      "Part. iHg",
                      "Lowest point"),
           col = c("black", "black", "black"),
           pch = c(16, 1, 17),
           cex = 0.9,
           bty = "n")
    
    
    
  } else {
    empty.plot("Missing MeHg data")
  }
  
  #### 8th plot: Plot out hgcA coverage ####
  if (RM.of.interest %in% MG.depth$RM) {
    
    MG.depth.info <- MG.depth %>%
      filter(RM == RM.of.interest)
    plot(x = MG.depth.info$hgcA_depth,
         y = MG.depth.info$depth,
         xlab = "hgcA depth of coverage", 
         ylab = "Depth (m)",
         xlim = c(0, 0.9),
         ylim = c(max.depth, 0),
         main = paste("Methylator abundance", 
                      sep = ""),
         col = "darkgreen",
         pch = 18)
  } else {
    empty.plot("No metagenomes here")
    
  }
  
  
}
dev.off()

