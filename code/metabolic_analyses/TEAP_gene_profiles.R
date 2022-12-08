#### code/metabolic_analyses/TEAP_gene_profiles.R ####
# Benjamin D. Peterson

# This script generates profiles of denitrification,
# sulfate reduction, and methanogenesis genes. 

#### Clean up crew on line 8 ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(ggpubr)
library(lubridate)
library(patchwork)
library(tidyverse)
source("code/HCC_plotting_needs.R")
rm(shape.vector)


#### Read in depth data ####
gene.data <- readRDS("dataEdited/metabolic_analyses/depth/metabolicProtein_depth_clean.rds") %>%
  filter(!is.na(geneName))


#### Save out summarized file ####
summarized.gene.data <- gene.data %>%
  group_by(date, RM, depth, elevation_m, geneName) %>%
  summarize(coverage = sum(coverage)) %>%
  spread(key = geneName,
         value = coverage)
write.csv(x = summarized.gene.data,
          file = "dataEdited/metabolic_analyses/summarized_gene_data.csv",
          row.names = FALSE)


#### Define function to plot multiple proteins ####
plot.profile.for.multiple.genes <- function(marker.depth.df,
                                            genesOfInterest,
                                            yearOfInterest,
                                            RMofInterest,
                                            gene.name.column = "geneName",
                                            xlab.to.use = "Gene coverage normalized\nto SCG coverage (100X)",
                                            depth_limits = c(80, 0),
                                            coverage_limits = NULL,
                                            show.mean.coverage = TRUE,
                                            color.vector.to.use = NULL,
                                            point.vector.to.use = NULL,
                                            line.vector.to.use = NULL,
                                            legend.position.to.use = "default",
                                            legend.title.to.use = element_blank(),
                                            titleToUse = NULL,
                                            DL = NULL,
                                            remove.depth.label = TRUE) {
  
  marker.depth.df[, "gene.names.to.use.for.filtering"] <- marker.depth.df[, gene.name.column]
  
  clean.coverage <- marker.depth.df %>%
    filter(gene.names.to.use.for.filtering %in% genesOfInterest) %>%
    filter(year(ymd(date)) == yearOfInterest) %>%
    filter(RM == RMofInterest) %>%
    group_by(gene.names.to.use.for.filtering, elevation_m) %>%
    summarise(coverage = sum(coverage)) %>%
    ungroup() %>%
    arrange(elevation_m)
  
  if (!is.null(DL)) {
    clean.coverage <- clean.coverage %>%
      spread(key = gene.names.to.use.for.filtering,
             value = coverage) %>%
      gather(key = gene.names.to.use.for.filtering,
             value = coverage,
             -elevation_m)
    clean.coverage[clean.coverage$coverage < 0.01, "coverage"] <- 0.010001
    
  }
  
  graph.of.interest <- clean.coverage %>%
    ggplot(aes(x = elevation_m,
               y = coverage)) +
    theme_classic() +
    scale_x_continuous(limits = depth_limits)+
    theme(axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))
  
  #### Add labels and title ####
  if (!is.null(titleToUse)) {
    graph.of.interest <- graph.of.interest +
      labs(title = titleToUse,
           y = xlab.to.use,
           x = "Elevation (m)")
    
  } else {
    graph.of.interest <- graph.of.interest +
      labs(y = xlab.to.use,
           x = "Elevation (m)")
  }
  
  #### Add colors if defined in the call ####
  if (!is.null(color.vector.to.use)) {
    graph.of.interest <- graph.of.interest +
      geom_line(aes(group = gene.names.to.use.for.filtering,
                    color = gene.names.to.use.for.filtering,
                    linetype = gene.names.to.use.for.filtering),
                size = 1) +
      geom_point(aes(group = gene.names.to.use.for.filtering,
                     color = gene.names.to.use.for.filtering,
                     shape = gene.names.to.use.for.filtering),
                 size = 2.5) +
      scale_colour_manual(values=color.vector.to.use)
  } else {
    graph.of.interest <- graph.of.interest +
      geom_point(shape = gene.names.to.use.for.filtering,
                 size = 2.5) +
      geom_line(aes(group = gene.names.to.use.for.filtering,
                    linetype = gene.names.to.use.for.filtering),
                size = 1)
  }
  
  #### Add point shapes if defined in the call ####
  if (!is.null(point.vector.to.use)) {
    graph.of.interest <- graph.of.interest +
      scale_shape_manual(values = point.vector.to.use)
  }
  #### Add line types if defined in the call ####
  if (!is.null(line.vector.to.use)) {
    graph.of.interest <- graph.of.interest +
      scale_linetype_manual(values = line.vector.to.use)
  }
  #### Constrain x-axis if defined in the call ####
  if (!is.null(coverage_limits)) {
    graph.of.interest <- graph.of.interest +
      scale_y_continuous(limits = coverage_limits)
  } else {
    max.coverage <- clean.coverage %>%
      select(coverage) %>%
      unlist() %>%
      max()
    graph.of.interest <- graph.of.interest +
      scale_y_continuous(limits = c(0.01, max.coverage))
    
  }
  
  if (show.mean.coverage == TRUE) {
    graph.of.interest <- graph.of.interest +
      stat_summary(geom = "point", fun = "mean",
                   col = "black", fill = "red",
                   size = 3, shape = 24)
  }
  
  #### Set up legend position ####
  if (legend.position.to.use[1] != "default") {
    graph.of.interest <- graph.of.interest +
      theme(legend.position = legend.position.to.use,
            legend.title = element_blank(),
            legend.key = element_rect(fill = "transparent", colour = "black"),
            legend.key.size = unit(1.75, 'lines'),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"))
  }
  
  if (remove.depth.label == TRUE) {
    graph.of.interest <- graph.of.interest +
      theme(axis.title.y = element_blank())
  }
  
  graph.of.interest + coord_flip()
  
}





#### TEAP profiles ####
color.vector.TEAP <- c(cb.translator["reddishpurple"],
                       cb.translator["blue"],
                       cb.translator["vermillion"])
names(color.vector.TEAP) <- c("narG", "dsrA", "mcrA")
shape.vector <- c(16, 17, 18)
names(shape.vector) <- c("narG", "dsrA", "mcrA")
line.vector <- c(1, 2, 3)
names(line.vector) <- c("narG", "dsrA", "mcrA")
N.max.coverage <- 60
TEAP.286.2017 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                                 genesOfInterest = names(color.vector.TEAP),
                                                 yearOfInterest = "2017",
                                                 RMofInterest = "286",
                                                 show.mean.coverage = FALSE,
                                                 remove.depth.label = FALSE,
                                                 depth_limits = c(545, 632),
                                                 coverage_limits = c(0, N.max.coverage),
                                                 color.vector.to.use = color.vector.TEAP,
                                                 point.vector.to.use = shape.vector,
                                                 line.vector.to.use = line.vector,
                                                 xlab.to.use = "Gene abundance (%)",
                                                 titleToUse = "2017 - RM286",
                                                 legend.position.to.use = c(0.6, 0.8))
TEAP.300.2017 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                                 genesOfInterest = names(color.vector.TEAP),
                                                 yearOfInterest = "2017",
                                                 RMofInterest = "300",
                                                 show.mean.coverage = FALSE,
                                                 remove.depth.label = FALSE,
                                                 depth_limits = c(545, 632),
                                                 coverage_limits = c(0, N.max.coverage),
                                                 point.vector.to.use = shape.vector,
                                                 color.vector.to.use = color.vector.TEAP,
                                                 line.vector.to.use = line.vector,
                                                 xlab.to.use = "Gene abundance (%)",
                                                 titleToUse = "2017 - RM300",
                                                 legend.position.to.use = c(0.6, 0.8))
TEAP.286.2018 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                                 genesOfInterest = names(color.vector.TEAP),
                                                 yearOfInterest = "2018",
                                                 RMofInterest = "286",
                                                 show.mean.coverage = FALSE,
                                                 remove.depth.label = FALSE,
                                                 depth_limits = c(545, 632),
                                                 coverage_limits = c(0, N.max.coverage),
                                                 color.vector.to.use = color.vector.TEAP,
                                                 point.vector.to.use = shape.vector,
                                                 line.vector.to.use = line.vector,
                                                 xlab.to.use = "Gene abundance (%)",
                                                 titleToUse = "2018 - RM286",
                                                 legend.position.to.use = c(0.6, 0.8))
TEAP.300.2018 <- plot.profile.for.multiple.genes(marker.depth.df = gene.data,
                                                 genesOfInterest = names(color.vector.TEAP),
                                                 yearOfInterest = "2018",
                                                 RMofInterest = "300",
                                                 show.mean.coverage = FALSE,
                                                 remove.depth.label = FALSE,
                                                 depth_limits = c(545, 632),
                                                 coverage_limits = c(0, N.max.coverage),
                                                 point.vector.to.use = shape.vector,
                                                 color.vector.to.use = color.vector.TEAP,
                                                 line.vector.to.use = line.vector,
                                                 xlab.to.use = "Gene abundance (%)",
                                                 titleToUse = "2018 - RM300",
                                                 legend.position.to.use = c(0.6, 0.8))



#### Arrange plots ####
pdf("results/metabolic_analyses/TEAP_gene_profiles.pdf",
    height = 7,
    width = 5)
ggarrange(TEAP.286.2017, TEAP.300.2017,
          TEAP.286.2018, TEAP.300.2018,
          nrow = 2,
          ncol = 2)
dev.off()


#### Generate separate plots to use in defining redox status of each sample ####
source("code/HCC_plotting_needs.R")
geochem.data <- readRDS("dataEdited/geochem/geochem_WC_adjusted_for_redox_classification.rds")

redox.status <- function(date.to.use,
                         RM.to.use) {
  geochem.data %>%
    filter(RM == RM.to.use,
           date == date.to.use) %>%
    ggplot(aes(x = 2,
               y = depth,
               color = redox_status)) +
    geom_point() +
    scale_color_manual(values = color.vector) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(title = "filler") +
    ylim(c(75, 0))
}
pdf("results/metabolic_analyses/TEAP_gene_profiles_REDOX_STATUS_TEMPLATE.pdf",
    height = 7,
    width = 5)
ggarrange(redox.status(date.to.use = "2017-09-25", 286), redox.status("2017-09-26", 300),
          redox.status("2018-09-24", 286), redox.status("2018-09-25", 300),
          nrow = 2,
          ncol = 2)
dev.off()
