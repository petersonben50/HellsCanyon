#### code/geochem/profile_functions.R #####



#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
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


# geochem.data.to.use = geochem.data
# date.range <- c("2018-04-10", "2018-04-13")
# RMs.to.use = NULL
# sulfide.plotting.factor = 100
# nitrate.plotting.factor = 5
# max.depth = 80


#### Function to plot redox ####
redox.plot <- function(geochem.data.to.use = geochem.data,
                       date.range,
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
                       max.depth = 80) {
  
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
    
    if (dim(geochem.data.site)[1] == 0) {
      empty.plot()
    } else {
      
      #### Add sulfide plot ####
      if(plot.Mn.instead.of.sulfide.YES.or.NO == "NO") {
        plot(x = geochem.data.site$f_inorganic_sulfide_mg_per_l*sulfide.plotting.factor,
             y = geochem.data.site$depth,
             xlab = "Sulfide (mg/L)",
             ylab = paste(year(date.range[1]), "-", month(date.range[1]), ":",
                          "Depth (m)",
                          sep = ""),
             xlim = c(0, 10),
             ylim = c(max.depth, 0),
             xaxt = "n",
             col = cb.translator[sulfide.color],
             pch = 18,
             cex = 1.2,
             main = paste("RM", RM.of.interest, sep = ""))
        lines(x = geochem.data.site$f_inorganic_sulfide_mg_per_l*sulfide.plotting.factor,
              y = geochem.data.site$depth,
              col = cb.translator[sulfide.color])
        # Add axis for sulfide
        axis(1,
             at = seq(0, 10, by = 2),
             labels = seq(0, 10, by = 2)/sulfide.plotting.factor)
        if (plot.Mn.with.sulfide == "YES") {
          points(x = geochem.data.site$f_mn_mg_per_l*Mn.plotting.factor,
                 y = geochem.data.site$depth,
                 col = cb.translator[Mn.color],
                 pch = 18)
          lines(x = geochem.data.site$f_mn_mg_per_l*Mn.plotting.factor,
                y = geochem.data.site$depth,
                col = cb.translator[Mn.color])
        }
        
      } else if (plot.Mn.instead.of.sulfide.YES.or.NO == "YES") {
        #### Add Mn plot instead of sulfide (if requested) ####
        plot(x = geochem.data.site$f_mn_mg_per_l*Mn.plotting.factor,
             y = geochem.data.site$depth,
             xlab = "Mn (mg/L)",
             ylab = paste(year(date.range[1]), "-", month(date.range[1]), ":",
                          "Depth (m)",
                          sep = ""),
             xlim = c(0, 10),
             ylim = c(max.depth, 0),
             xaxt = "n",
             col = cb.translator[Mn.color],
             pch = 18,
             cex = 1.2,
             main = paste("RM", RM.of.interest, sep = ""))
        lines(x = geochem.data.site$f_mn_mg_per_l*Mn.plotting.factor,
              y = geochem.data.site$depth,
              col = cb.translator[Mn.color])
        # Add axis for sulfide
        axis(1,
             at = seq(0, 10, by = 2),
             labels = seq(0, 10, by = 2)/Mn.plotting.factor)
        
        
      }
      
      if (use.nitrate.chloride.ratios.YES.or.NO == "YES") {
        points(x = geochem.data.site$f_no3_mg_n_per_l/geochem.data.site$f_cl_mg_per_l*nitrate.chloride.plotting.factor,
               y = geochem.data.site$depth,
               col = cb.translator[nitrate.color],
               pch = 18)
        lines(x = geochem.data.site$f_no3_mg_n_per_l/geochem.data.site$f_cl_mg_per_l*nitrate.chloride.plotting.factor,
              y = geochem.data.site$depth,
              col = cb.translator[nitrate.color])
        # Add label for nitrate
        title(xlab = "Nitrate/Chloride",
              line = 4.5)
      } else if (use.nitrate.chloride.ratios.YES.or.NO == "NO") {
        points(x = geochem.data.site$f_no3_mg_n_per_l*nitrate.plotting.factor,
               y = geochem.data.site$depth,
               col = cb.translator[nitrate.color],
               pch = 18)
        lines(x = geochem.data.site$f_no3_mg_n_per_l*nitrate.plotting.factor,
              y = geochem.data.site$depth,
              col = cb.translator[nitrate.color])
        # Add label for nitrate
        title(xlab = "Nitrate (mgN/L)",
              line = 4.5)
        
      }
      
      points(x = geochem.data.site$MeHg_diss_ngL*MeHg.plotting.factor,
             y = geochem.data.site$depth,
             col = cb.translator[MeHg.color],
             pch = 18)
      lines(x = geochem.data.site$MeHg_diss_ngL*MeHg.plotting.factor,
            y = geochem.data.site$depth,
            col = cb.translator[MeHg.color])

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


