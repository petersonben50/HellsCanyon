#### code/geochem/profile_functions.R #####



#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(viridisLite)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")



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

#### Test data ####
# date.list <- list(
#   c("2015-07-14", "2015-07-15"),
#   c("2015-08-10", "2015-08-11"),
#   c("2015-09-08", "2015-09-09"),
#   c("2015-10-19", "2015-10-19")
# )
# 
# geochem.data.to.use = geochem.data
# DO.data.to.use = DO.data
# date.range
# DO.date = NULL
# RMs.to.use = NULL
# plot.Mn.instead.of.sulfide.YES.or.NO = "NO"
# sulfide.plotting.factor = 100
# sulfide.color = "blue"
# plot.Mn.with.sulfide = "YES"
# Mn.plotting.factor = 3
# Mn.color = "orange"
# use.nitrate.chloride.ratios.YES.or.NO = "NO"
# nitrate.plotting.factor = 5
# nitrate.chloride.plotting.factor = 100
# nitrate.color = "reddishpurple"
# MeHg.plotting.factor = 3
# MeHg.color = "bluishgreen"
# DO.plotting.factor = 4
# DO.color = "black"
# min.elevation = 545
# max.elevation = 632
# plot.DO = "NO"
# geochem.data.to.use = geochem.data
# date.range = c("2017-09-25", "2017-09-28")
# RMs.to.use = c(286, 300, 310, 318)
# nitrate.plotting.factor = nitrate.plotting.factor.to.use
# MeHg.plotting.factor = 2.5
# plot.Mn.instead.of.sulfide.YES.or.NO = "NO"
# Mn.plotting.factor = Mn.plotting.factor.to.use
# plot.DO = "YES"
# 

#### Function to plot redox ####
redox.plot <- function(geochem.data.to.use = geochem.data,
                       DO.data.to.use = DO.data,
                       date.range,
                       DO.date = NULL,
                       RMs.to.use = NULL,
                       plot.Mn.instead.of.sulfide.YES.or.NO = "NO",
                       sulfide.plotting.factor = 100,
                       sulfide.color = "blue",
                       plot.Mn.with.sulfide = "YES",
                       Mn.plotting.factor = 3,
                       Mn.color = "orange",
                       use.nitrate.chloride.ratios.YES.or.NO = "NO",
                       nitrate.plotting.factor = 5,
                       nitrate.chloride.plotting.factor = 100,
                       nitrate.color = "reddishpurple",
                       MeHg.plotting.factor = 3,
                       MeHg.color = "bluishgreen",
                       DO.plotting.factor = 4,
                       DO.color = "black",
                       min.elevation = 545,
                       max.elevation = 632,
                       plot.DO = "NO") {
  
  geochem.data.date <- geochem.data.to.use %>%
    spread(key = constituent,
           value = concentration) %>%
    filter(date >= as.Date(date.range[1]) &
             date <= as.Date(date.range[2]))
  
  if (is.null(RMs.to.use)){
    RMs.to.use <- sort(unique(geochem.data.date$RM))
    }

  for (RM.of.interest in RMs.to.use) {
    
    geochem.data.site <- geochem.data.date %>%
      filter(RM == RM.of.interest) %>%
      arrange(depth)
    date.of.sampling.at.RM <- geochem.data.site %>%
      select(date) %>%
      unlist(use.names = FALSE) %>%
      unique()
    
    if (dim(geochem.data.site)[1] == 0) {
      empty.plot()
    } else {
      
      
      
      #### Add sulfide plot ####
      if(plot.Mn.instead.of.sulfide.YES.or.NO == "NO") {
        plot(x = geochem.data.site$f_inorganic_sulfide_mg_per_l*sulfide.plotting.factor,
             y = geochem.data.site$elevation_m,
             xlab = "Sulfide (mg/L)",
             ylab = paste(date.of.sampling.at.RM, "\n",
                          "Depth (m)",
                          sep = ""),
             xlim = c(0, 10),
             ylim = c(min.elevation, max.elevation),
             xaxt = "n",
             col = cb.translator[sulfide.color],
             pch = 18,
             cex = 1.2,
             main = paste("RM", RM.of.interest, sep = ""))
        lines(x = geochem.data.site$f_inorganic_sulfide_mg_per_l*sulfide.plotting.factor,
              y = geochem.data.site$elevation_m,
              col = cb.translator[sulfide.color])
        # Add axis for sulfide
        axis(1,
             at = seq(0, 10, by = 2),
             labels = seq(0, 10, by = 2)/sulfide.plotting.factor)
        if (plot.Mn.with.sulfide == "YES") {
          points(x = geochem.data.site$f_mn_mg_per_l*Mn.plotting.factor,
                 y = geochem.data.site$elevation_m,
                 col = cb.translator[Mn.color],
                 pch = 18)
          lines(x = geochem.data.site$f_mn_mg_per_l*Mn.plotting.factor,
                y = geochem.data.site$elevation_m,
                col = cb.translator[Mn.color])
        }
        
      } else if (plot.Mn.instead.of.sulfide.YES.or.NO == "YES") {
        #### Add Mn plot instead of sulfide (if requested) ####
        plot(x = geochem.data.site$f_mn_mg_per_l*Mn.plotting.factor,
             y = geochem.data.site$elevation_m,
             xlab = "Mn (mg/L)",
             ylab = paste(date.of.sampling.at.RM, "\n",
                          "Depth (m)",
                          sep = ""),
             xlim = c(0, 10),
             ylim = c(min.elevation, max.elevation),
             xaxt = "n",
             col = cb.translator[Mn.color],
             pch = 18,
             cex = 1.2,
             main = paste("RM", RM.of.interest, sep = ""))
        lines(x = geochem.data.site$f_mn_mg_per_l*Mn.plotting.factor,
              y = geochem.data.site$elevation_m,
              col = cb.translator[Mn.color])
        # Add axis for Mn
        axis(1,
             at = seq(0, 10, by = 2),
             labels = seq(0, 10, by = 2)/Mn.plotting.factor)
        
        
      }
      
      if (use.nitrate.chloride.ratios.YES.or.NO == "YES") {
        points(x = geochem.data.site$f_no3_mg_n_per_l/geochem.data.site$f_cl_mg_per_l*nitrate.chloride.plotting.factor,
               y = geochem.data.site$elevation_m,
               col = cb.translator[nitrate.color],
               pch = 18)
        lines(x = geochem.data.site$f_no3_mg_n_per_l/geochem.data.site$f_cl_mg_per_l*nitrate.chloride.plotting.factor,
              y = geochem.data.site$elevation_m,
              col = cb.translator[nitrate.color])
        # Add label for nitrate/chloride ratio
        title(xlab = "Nitrate/Chloride",
              line = 4.5)
      } else if (use.nitrate.chloride.ratios.YES.or.NO == "NO") {
        points(x = geochem.data.site$f_no3_mg_n_per_l*nitrate.plotting.factor,
               y = geochem.data.site$elevation_m,
               col = cb.translator[nitrate.color],
               pch = 18)
        lines(x = geochem.data.site$f_no3_mg_n_per_l*nitrate.plotting.factor,
              y = geochem.data.site$elevation_m,
              col = cb.translator[nitrate.color])
        # Add label for nitrate
        title(xlab = "Nitrate (mgN/L)",
              line = 4.5)
        
      }
      
      points(x = geochem.data.site$MeHg_diss_ngL*MeHg.plotting.factor,
             y = geochem.data.site$elevation_m,
             col = cb.translator[MeHg.color],
             pch = 18)
      lines(x = geochem.data.site$MeHg_diss_ngL*MeHg.plotting.factor,
            y = geochem.data.site$elevation_m,
            col = cb.translator[MeHg.color])

      
      
      if (plot.DO == "YES") {
        
        if (is.null(DO.date)) {
          DO.data.date <- DO.data.to.use %>%
            filter(date == date.of.sampling.at.RM &
                     RM == RM.of.interest)
        } else {
          DO.data.date <- DO.data.to.use %>%
            filter(date == DO.date,
                   RM == RM.of.interest)
        }
        
        points(x = DO.data.date$diss_oxy_mg_per_l*0.5,
               y = DO.data.date$elevation_m,
               col = cb.translator[DO.color],
               pch = 18)
      }
      # Add axis for nitrate
      axis(1,
           line = 3,
           at = seq(0, 10, by = 2),
           labels = seq(0, 10, by = 2)/nitrate.plotting.factor)
      # Add axis for MeHg
      axis(1,
           line = 6,
           at = seq(0, 10, by = 2),
           labels = seq(0, 10, by = 2)/MeHg.plotting.factor)

      # Add label for MeHg
      title(xlab = "MeHg (ng/L)",
            line = 7.5)
      }
  }
}




#### Function to plot time course over year of single parameter ####
time.course.profile.plot <- function(geochem.data.to.use,
                                     years.to.use,
                                     RMs.to.use,
                                     parameter.to.plot,
                                     xlabel.to.use,
                                     color.ramp.to.use = NULL,
                                     color.vector.to.use = NULL,
                                     concentrations.to.use = c(0, 2.2)
                                     ) {
  geochem.data.to.use.temporary <- geochem.data.to.use %>%
    filter(year(date) %in% years.to.use,
           RM %in% RMs.to.use,
           constituent == parameter.to.plot) %>%
    arrange(date, elevation_m)
  
  if (is.null(color.vector.to.use)) {
    if (is.null(color.ramp.to.use)) {
      color.ramp.to.use = c("blue", "red")
    }
    colorize.function = colorRampPalette(color.ramp.to.use)
    color.vector.to.use = colorize.function(length(unique(month(geochem.data.to.use.temporary$date))))
  }
  geochem.data.to.use.temporary %>%
    ggplot(aes(x = concentration,
               y = elevation_m,
               col = month(date,label = TRUE),
               group = date)) +
    geom_point() +
    geom_path() +
    theme_classic() +
    scale_color_manual(values = color.vector.to.use,
                       name = "Month") +
    ylim(c(543, 635)) +
    xlim(concentrations.to.use) +
    labs(x = xlabel.to.use,
         y = "Elevation (m)",
         title = years.to.use)
    
}




#### Function to plot time course over year of Seabird parameter ####
seabird.time.course.profile.plot <- function(seabird.data.to.use,
                                             years.to.use,
                                             RMs.to.use,
                                             parameter.to.plot,
                                             xlabel.to.use,
                                             color.ramp.to.use = NULL,
                                             concentrations.to.use = c(0, 2.2)
) {
  seabird.data.to.use.temporary <- seabird.data.to.use %>%
    gather(key = constituent,
           value = concentration,
           -c(1:4)) %>%
    filter(year(date) %in% years.to.use,
           RM %in% RMs.to.use,
           constituent == parameter.to.plot) %>%
    arrange(date, elevation_m) %>%
    mutate(concentration = as.numeric(concentration))
  
  if (is.null(color.ramp.to.use)) {
    color.ramp.to.use <- c("blue", "red")
  }
  colorize.function <- colorRampPalette(color.ramp.to.use)
  
  seabird.data.to.use.temporary %>%
    ggplot(aes(x = concentration,
               y = elevation_m,
               col = month(date,label = TRUE),
               group = date)) +
    geom_point() +
    geom_path() +
    theme_classic() +
    scale_color_manual(values = colorize.function(length(unique(month(seabird.data.to.use.temporary$date,
                                                                      label = TRUE)))),
                       name = "Month") +
    ylim(c(543, 635)) +
    xlim(concentrations.to.use) +
    labs(x = xlabel.to.use,
         y = "Elevation (m)",
         title = years.to.use)
  
}

