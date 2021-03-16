#### code/profiles/2018_intensive.R ####
# Benjamin D. Peterson

#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(naniar)
library(readxl)
library(tidyverse)


#### Read in geochem data ####
geochem.data <- read.csv("dataEdited/waterChemistry/geochem_WC_2015_2018.csv",
                         stringsAsFactors = FALSE) %>%
  filter(year(date) == "2018") %>%
  filter(month(date) == "9") %>%
  filter(day(date) >= 23) %>%
  mutate(depth = as.numeric(depth)) %>%
  arrange(date, RM, depth)
# Replace non-detects with zeroes in sulfide column
geochem.data$sulfide_mg.L[geochem.data$sulfide_mg.L == "nd"] <- 0
geochem.data$sulfide_mg.L <- as.numeric(geochem.data$sulfide_mg.L)


#### Read in Hg data ####
MeHg.data <- read.csv("dataEdited/waterChemistry/Hg_2015_2018.csv",
                      stringsAsFactors = FALSE) %>%
  filter(year(date) == "2018") %>%
  filter(month(date) == "9") %>%
  filter(day(date) >= 23) %>%
  mutate(depth = as.numeric(depth)) %>%
  arrange(date, RM, depth)


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
  filter(year(date) == 2018)


#### List of sites ####
site.list <- sort(unique(c(geochem.data$RM,
                           MeHg.data$RM)))


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



#### Plot temp/DO, redox, and MeHg in separate plots data ####

pdf("results/geochem/profiles/2018_summary_profiles.pdf",
    height = 15,
    width = 15)

par(mfrow = c(4, 8),
    mgp=c(1.5,0.4,0),
    tck=-0.008)

for (RM.of.interest in site.list) {


  par(mar = c(6, 4.5, 3, 1))
  
  #### First plot: Anions ####
  
  if (RM.of.interest %in% geochem.data$RM) {
    
    geochem.data.site <- geochem.data %>%
      filter(RM == RM.of.interest) %>%
      replace_with_na_all(condition = ~.x == "na") %>%
      select(c(depth, SO4_mg.L, NO3_mgN.L)) %>%
      as.data.frame()
    
    max.depth <- max(geochem.data.site$depth)
    
    plot(x = geochem.data.site$SO4_mg.L/5,
         y = geochem.data.site$depth,
         xlab = "Sulfate (mg/L)",
         ylab = "Depth (m)",
         xlim = c(0, 10),
         ylim = c(max.depth, 0),
         xaxt = "n",
         col = "blue",
         pch = 18,
         cex = 1.2,
         main = "Anion concentrations")
    points(x = geochem.data.site$NO3_mgN.L*4,
           y = geochem.data.site$depth,
           col = "red",
           pch = 18)
    
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
    
    
    # Add DNA sampling information
    
    NA.filters.site <- MG.depth %>%
      filter(RM == RM.of.interest)
    
    if (length(NA.filters.site$depth) > 0) {
      points(x = rep(10, length(NA.filters.site$depth)),
             y = NA.filters.site$depth,
             pch = 18,
             col = "green")
      text(x = rep(7, length(NA.filters.site$depth)),
           y = NA.filters.site$depth,
           labels = NA.filters.site$met,
           col = "green")
    }
    
    
    
    

    # Add legend
    legend("topleft",
           legend = c("Sulfate (mg/L)",
                      "Nitrate (mgN/L)",
                      "Nitrite (ugN/L)"),
           col = c("blue",
                   "red",
                   "red"),
           pch = c(18, 18, 5),
           bty = "n")
    
    # Clean up
    rm(geochem.data.site)
    
  } else {
    empty.plot("No anion data")
  }
  
  title(ylab = paste("RM",
                     RM.of.interest,
                     sep = ""),
        line = 3,
        cex.lab = 2)
  
  par(mar = c(6, 3, 3, 1))
  
  
  #### Second plot: Anion ratios ####
  if (RM.of.interest %in% geochem.data$RM) {
    
    geochem.data.site <- geochem.data %>%
      filter(RM == RM.of.interest) %>%
      replace_with_na_all(condition = ~.x == "na") %>%
      mutate(SO4.Cl.ratio = SO4_mg.L / Cl_mg.L) %>%
      mutate(NO3.Cl.ratio = NO3_mgN.L / Cl_mg.L) %>%
      select(c(depth, SO4_mg.L, NO3_mgN.L,
               SO4.Cl.ratio, NO3.Cl.ratio)) %>%
      as.data.frame()
    
    plot(x = geochem.data.site$SO4.Cl.ratio,
         y = geochem.data.site$depth,
         xlab = "Sulfate/Chloride",
         ylab = "Depth (m)",
         xlim = c(0, 10),
         ylim = c(max.depth, 0),
         col = "blue",
         pch = 18,
         cex = 1.2,
         main = "Anion ratios")
    points(x = geochem.data.site$NO3.Cl.ratio*100,
           y = geochem.data.site$depth,
           col = "red",
           pch = 18)

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
           legend = c("Sulfate/Chloride",
                      "Nitrate/Chloride"),
           col = c("blue",
                   "red"),
           pch = c(18, 18, 5),
           bty = "n")
    
    # Clean up
    rm(geochem.data.site)   
    
    } else {
    empty.plot("No anion data")
    }
  
  
  #### Third and fourth plot: Metal data ####
  if (RM.of.interest %in% geochem.data$RM) {
    
    geochem.data.site <- geochem.data %>%
      filter(RM == RM.of.interest) %>%
      replace_with_na_all(condition = ~.x == "na") %>%
      select(c(depth, Fe_ug.L, Mn_ug.L, Fe.part_ug.L,
               Mn.part_ug.L)) %>%
      as.data.frame() %>%
      arrange(depth)
    
    iron.max <- max(c(geochem.data$Fe_ug.L,
                      geochem.data$Fe.part_ug.L - geochem.data$Fe_ug.L),
                    na.rm = TRUE)
    
    plot(x = geochem.data.site$Fe_ug.L,
         y = geochem.data.site$depth,
         xlab = "Fe (ug/L)",
         ylab = "Depth (m)",
         xlim = c(0, iron.max),
         ylim = c(max.depth, 0),
         pch = 18,
         cex = 0,
         main = "Iron concentrations")

    # Add in dissolved Fe
    points(x = geochem.data.site$Fe_ug.L,
           y = geochem.data.site$depth,
           col = "orange",
           pch = 18)
    # Add in particulate Fe
    points(x = geochem.data.site$Fe.part_ug.L - geochem.data.site$Fe_ug.L,
           y = geochem.data.site$depth,
           col = "orange",
           pch = 5)
    
    # Add legend
    legend("topright",
           legend = c("Fe (diss.)",
                      "Fe (part.)"),
           col = c("orange",
                   "orange"),
           pch = c(18, 5),
           bty = "n")
    
    
    
    Mn.max <- max(c(geochem.data$Mn_ug.L,
                      geochem.data$Mn.part_ug.L - geochem.data$Mn_ug.L),
                    na.rm = TRUE)
    
    plot(x = geochem.data.site$Mn_ug.L,
         y = geochem.data.site$depth,
         xlab = "Mn (ug/L)",
         ylab = "Depth (m)",
         xlim = c(0, Mn.max),
         ylim = c(max.depth, 0),
         pch = 18,
         cex = 0,
         main = "Mn concentrations")
    
    # Add in particulate Mn
    points(x = geochem.data.site$Mn.part_ug.L - geochem.data.site$Mn_ug.L,
           y = geochem.data.site$depth,
           col = "purple",
           pch = 5)
    lines(x = geochem.data.site$Mn.part_ug.L - geochem.data.site$Mn_ug.L,
          y = geochem.data.site$depth,
          col = "purple",
          lwd = 0.8,
          lty = 2)
    
    # Add in dissolved Mn
    points(x = geochem.data.site$Mn_ug.L,
           y = geochem.data.site$depth,
           col = "purple",
           pch = 18)
    lines(x = geochem.data.site$Mn_ug.L,
          y = geochem.data.site$depth,
          col = "purple",
          lwd = 0.8)
    
    # Add legend
    legend("topright",
           legend = c("Mn (diss.)",
                      "Mn (part.)"),
           col = c("purple",
                   "purple"),
           pch = c(18, 5),
           bty = "n")
    
    # Clean up
    rm(geochem.data.site)   
    
  } else {
    empty.plot("No anion data")
  }
  
  
  #### Fifth plot: sulfide data ####
  if (RM.of.interest %in% geochem.data$RM) {
    
    geochem.data.site <- geochem.data %>%
      filter(RM == RM.of.interest) %>%
      replace_with_na_all(condition = ~.x == "na") %>%
      select(c(depth, sulfide_mg.L)) %>%
      as.data.frame() %>%
      arrange(depth)
    
    sulfide.max <- max(geochem.data.site$sulfide_mg.L,
                       na.rm = TRUE)
    if (sulfide.max <= 1) {
      sulfide.max <- 1
    }
    
    plot(x = geochem.data.site$sulfide_mg.L,
         y = geochem.data.site$depth,
         xlab = "Sulfide (mg/L)",
         ylab = "Depth (m)",
         xlim = c(0, sulfide.max),
         ylim = c(max.depth, 0),
         pch = 18,
         cex = 0,
         main = "Sulfide concentrations")
    
    # Add in sulfide points
    points(x = geochem.data.site$sulfide_mg.L,
           y = geochem.data.site$depth,
           col = "blue",
           pch = 18)
    lines(x = geochem.data.site$sulfide_mg.L,
          y = geochem.data.site$depth,
          col = "blue",
          lwd = 0.8,
          lty = 1)
    
    # Add legend
    legend("topright",
           legend = c("Sulfide"),
           col = c("blue",
                   "orange"),
           pch = c(18),
           bty = "n")
    # Clean up
    rm(geochem.data.site)   
    
  } else {
    empty.plot("No sulfide data")
  }
  

  #### Fifth and sixth plot: MeHg ####
  if (RM.of.interest %in% MeHg.data$RM) {
    
    MeHg.data.site <- MeHg.data %>%
      filter(RM == RM.of.interest)

    plot(x = MeHg.data.site$FMHG,
         y = MeHg.data.site$depth,
         xlab = "MeHg (ng/L)", 
         ylab = "Depth (m)",
         xlim = c(0, 2.5),
         ylim = c(max.depth, 0),
         main = paste("Hg data", 
                      sep = ""),
         col = "firebrick3",
         pch = 18)
    lines(x = MeHg.data.site$FMHG,
          y = MeHg.data.site$depth,
          col = "firebrick3",
          lwd = 0.8)
    
    # Add in particulate MeHg
    points(x = MeHg.data.site$PMHG,
           y = MeHg.data.site$depth,
           col = "firebrick3",
           pch = 5)
    # Add in particulate MeHg
    lines(x = MeHg.data.site$PMHG,
          y = MeHg.data.site$depth,
          col = "firebrick3",
          lwd = 0.8,
          lty = 2)
    
    legend("topright",
           legend = c("Diss. MeHg",
                      "Part. MeHg"),
           col = c("firebrick3",
                   "firebrick3"),
           pch = c(18,
                   5),
           bty = "n") 
    
    
    plot(x = MeHg.data.site$FTHG - MeHg.data.site$FMHG,
         y = MeHg.data.site$depth,
         xlab = "Hg (ng/L)", 
         ylab = "Depth (m)",
         xlim = c(0, 2),
         ylim = c(max.depth, 0),
         main = paste("iHg data", 
                      sep = ""),
         col = "black",
         pch = 18)
    lines(x = MeHg.data.site$FTHG - MeHg.data.site$FMHG,
          y = MeHg.data.site$depth,
          col = "black",
          pch = 5,
          lwd = 0.8)

    # Add in particulate THg
    points(x = MeHg.data.site$PTHG - MeHg.data.site$PMHG,
           y = MeHg.data.site$depth,
           col = "black",
           pch = 5)
    lines(x = MeHg.data.site$PTHG - MeHg.data.site$PMHG,
          y = MeHg.data.site$depth,
          col = "black",
          pch = 5,
          lwd = 0.8,
          lty = 2)
    
    legend("topright",
           legend = c("Diss. iHg",
                      "Part. iHg"),
           col = c("black",
                   "black"),
           pch = c(18,
                   5),
           bty = "n")
  } else {
    empty.plot("Missing MeHg data")
  }
  
  #### 7th plot: Plot out hgcA coverage ####
  if (RM.of.interest %in% MG.depth$RM) {
    
    MG.depth.info <- MG.depth %>%
      filter(RM == RM.of.interest)
    plot(x = MG.depth.info$hgcA_depth,
         y = MG.depth.info$depth,
         xlab = "hgcA depth of coverage", 
         ylab = "Depth (m)",
         xlim = c(0, 3.3),
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

