#### code/gene_plotting_functions.R ####
# Written by Benjamin D. Peterson



#### Sonde profile function ####
sonde.profile <- function(seabird.data.of.interest,
                          RM.of.interest,
                          year.of.interest,
                          depth.range.of.interest,
                          concentration.of.interest,
                          legend.location.of.interest,
                          scaling.vector.to.use = NULL,
                          color.vector.to.use) {
  
  if (!is.null(scaling.vector.to.use)) {
    for (modification.position in 1:length(scaling.vector.to.use)) {
      modifications.to.keep <- seabird.data.of.interest %>%
        filter(constituent == names(scaling.vector.to.use)[modification.position]) %>%
        mutate(concentration = concentration * scaling.vector.to.use[modification.position])
      seabird.data.of.interest <- seabird.data.of.interest %>%
        filter(constituent != names(scaling.vector.to.use)[modification.position]) %>%
        rbind(modifications.to.keep)
    }
  }

  plot.to.plot <- seabird.data.of.interest %>%
    filter(RM == RM.of.interest,
           year(date) == year.of.interest,
           constituent %in% names(color.vector.to.use)) %>%
    ggplot(aes(x = depth,
               y = concentration,
               color = constituent)) +
    geom_point(aes(color = constituent),
               size = 0.5) +
    geom_line(aes(color = constituent,
                  linetype = constituent),
              size = 0.5) +
    scale_colour_manual(values = color.vector.to.use,
                        labels = labels.vector) +
    scale_linetype_manual(values = line.vector,
                          labels = labels.vector) +
    coord_flip(xlim = depth.range.of.interest,
               ylim = concentration.of.interest) +
    ylab("Constituent") +
    xlab("Depth (m)") +
    theme_classic() +
    theme(legend.position = legend.location.of.interest,
          legend.title = element_blank(),
          legend.key = element_rect(fill = "transparent", colour = "black"),
          legend.key.size = unit(1.75, 'lines'),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black")) #+ 
    # scale_y_continuous(sec.axis = sec_axis(~ . * scaling.vector.to.use[1],
    #                                        name = labels.vector[names(scaling.vector.to.use)[1]]))
  
  plot.to.plot
  
}



#### Set plotting function to generate plot of geochem data ####
plot.redox.profile <- function(geochem.data.of.interest = geochem.data,
                               RM.of.interest,
                               date.of.interest,
                               depth.range.of.interest,
                               concentration.of.interest,
                               legend.location.of.interest,
                               remove.depth.label = TRUE,
                               scaling.vector.to.use = NULL) {
  geochem.data.of.interest <- geochem.data.of.interest %>%
    filter(RM == RM.of.interest,
           date == as.Date(date.of.interest)) %>%
    group_by(RM, depth, date, constituent) %>%
    summarize(concentration = mean(concentration)) %>%
    filter(constituent %in% names(color.vector))

  if (!is.null(scaling.vector.to.use)) {
    for (modification.position in 1:length(scaling.vector.to.use)) {
      modifications.to.keep <- geochem.data.of.interest %>%
        filter(constituent == names(scaling.vector.to.use)[modification.position]) %>%
        mutate(concentration = concentration * scaling.vector.to.use[modification.position])
      geochem.data.of.interest <- geochem.data.of.interest %>%
        filter(constituent != names(scaling.vector.to.use)[modification.position]) %>%
        rbind(modifications.to.keep)
    }
  }

  making.a.plot <- geochem.data.of.interest %>%
    ggplot(aes(x = depth,
               y = concentration,
               color = constituent,
               linetype = constituent,
               shape = constituent)) +
    geom_point(aes(color = constituent),
               size = 1.5) +
    geom_line(aes(color = constituent),
              size = 0.6) +
    scale_colour_manual(values = color.vector,
                        labels = labels.vector) +
    scale_shape_manual(values = points.vector,
                       labels = labels.vector) +
    scale_linetype_manual(values = line.vector,
                          labels = labels.vector) +
    coord_flip(xlim = depth.range.of.interest,
               ylim = concentration.of.interest) +
    ylab("Concentration (mg/L)") +
    xlab("Depth (m)") +
    theme_classic() +
    theme(legend.position = legend.location.of.interest,
          legend.title = element_blank(),
          legend.key = element_rect(fill = "transparent", colour = "black"),
          legend.key.size = unit(1.75, 'lines'),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"))
  
  if (remove.depth.label == TRUE) {
    making.a.plot <- making.a.plot +
      theme(axis.title.y = element_blank())
  }
  
  
  making.a.plot
}













#### Define function to plot hgcA and MeHg ####

plot.hgcA.and.MeHg <- function(data.to.use,
                               date.of.interest,
                               RM.of.interest,
                               xlab.to.use = "hgcA abundance (%), MeHg (ng/L)",
                               depth.range.of.interest,
                               legend.location.of.interest,
                               remove.depth.label = TRUE,
                               abundance.limits) {
  graph.of.interest <- data.to.use %>%
    filter(date == date.of.interest,
           RM == RM.of.interest) %>%
    ggplot(aes(x = depth,
               y = amount,
               color = constituent,
               linetype = constituent,
               shape = constituent)) +
    geom_point(aes(color = constituent),
               size = 2) +
    geom_line(aes(color = constituent),
              size = 1) +
    scale_colour_manual(values = color.vector,
                        labels = labels.vector) +
    scale_shape_manual(values = points.vector,
                       labels = labels.vector) +
    scale_linetype_manual(values = line.vector,
                          labels = labels.vector) +
    scale_y_continuous(trans = 'log10',
                       limits = abundance.limits) +
    coord_flip(xlim = depth.range.of.interest) +
    theme_classic()  +
    ylab(xlab.to.use) +
    theme(legend.position = legend.location.of.interest,
          legend.title = element_blank(),
          legend.key = element_rect(fill = "transparent", colour = "black"),
          legend.key.size = unit(1.75, 'lines'),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"))
  
  if (remove.depth.label == TRUE) {
    graph.of.interest <- graph.of.interest +
      theme(axis.title.y = element_blank())
  }
  
  graph.of.interest
  
  }



#### Define function to plot one gene with taxonomy ####
# data.to.use dataframe must have column with date, RM,
# depth of sample collection, coverage of gene, and
# a column with a taxonomic identifier.
plot.profile.of.gene.with.taxonomy <- function(data.to.use,
                                               RM.of.interest,
                                               year.of.interest,
                                               legend.position.to.use = "default",
                                               depth_limits = c(80, 0),
                                               coverage_limits = NULL,
                                               taxonomy.column.name,
                                               color.vector.to.use = color.vector,
                                               xlab.to.use = "Depth (m)",
                                               ylab.to.use = "Gene coverage (per 100X SCG coverage)",
                                               titleToUse = element_blank()) {
  data.to.use[, "taxonomy"] <- data.to.use[, taxonomy.column.name]

  graph.to.make <- data.to.use %>%
    filter(year(date) == year.of.interest,
           RM == RM.of.interest) %>%
    ggplot(aes(x = depth,
               y = coverage,
               fill = taxonomy)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color.vector.to.use) +
    xlim(depth_limits)+
    theme_classic()

  #### Constrain x-axis if defined in the call ####
  if (!is.null(coverage_limits)) {
    graph.to.make <- graph.to.make +
      scale_y_continuous(limits = coverage_limits)
  } else {
    max.coverage <- data.to.use %>%
      filter(year(date) == year.of.interest,
             RM == RM.of.interest) %>%
      group_by(depth) %>%
      summarise(coverage = sum(coverage)) %>%
      select(coverage) %>%
      unlist() %>%
      max()
    graph.to.make <- graph.to.make +
      scale_y_continuous(limits = c(0, max.coverage))

  }

  if (legend.position.to.use[1] != "default") {
    graph.to.make <- graph.to.make +
      theme(legend.position = legend.position.to.use,
            legend.title = element_blank())
  }

  #### Add labels and title ####
  graph.to.make <- graph.to.make +
    labs(x = xlab.to.use,
         y = ylab.to.use,
         title = titleToUse)

  graph.to.make <- graph.to.make +
    theme(axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"))

  graph.to.make +
    coord_flip()

}








#### Define function to plot one gene with taxonomy ####
# data.to.use dataframe must have column with date, RM,
# depth of sample collection, coverage of gene, and
# a column with a taxonomic identifier.
plot.profile.of.gene.with.taxonomy.function.noDepth <- function(data.to.use = hgcA.data,
                                                                RM.of.interest = 286,
                                                                year.of.interest = 2017,
                                                                legend.position.to.use = "default",
                                                                coverage_limits = NULL,
                                                                taxonomy.column.name = "manual_classification",
                                                                color.vector.tax.to.use = color.vector.tax,
                                                                ylab.to.use = "Gene coverage (per 100X SCG coverage)",
                                                                xlab.to.use = "Depth (m)",
                                                                titleToUse = element_blank()) {
  data.to.use[, "taxonomy"] <- data.to.use[, taxonomy.column.name]

  data.to.use <- data.to.use %>%
    filter(year(date) == year.of.interest,
           RM == RM.of.interest) %>%
    mutate(depth = as.character(depth))


  depth.order <- unique(data.to.use$depth)[order(as.numeric(unique(data.to.use$depth)),
                                                 decreasing = TRUE)]
  tax.order.by.function <- data.to.use %>%
    arrange(predicted_metabolism) %>%
    select(manual_classification) %>%
    unlist(use.names = FALSE) %>%
    unique() %>%
    as.character()

  data.to.use <- data.to.use %>%
    mutate(depth = fct_relevel(depth,
                               depth.order),
           taxonomy = fct_relevel(taxonomy,
                                  tax.order.by.function))
  graph.to.make <- data.to.use %>%
    ggplot(aes(x = depth,
               y = coverage,
               fill = taxonomy)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color.vector.tax.to.use) +
    theme_classic()

  #### Constrain x-axis if defined in the call ####
  if (!is.null(coverage_limits)) {
    graph.to.make <- graph.to.make +
      scale_y_continuous(limits = coverage_limits,
                         expand = c(0, 0))
  } else {
    max.coverage <- data.to.use %>%
      filter(year(date) == year.of.interest,
             RM == RM.of.interest) %>%
      group_by(depth) %>%
      summarise(coverage = sum(coverage)) %>%
      select(coverage) %>%
      unlist() %>%
      max()
    graph.to.make <- graph.to.make +
      scale_y_continuous(limits = c(0, max.coverage),
                         expand = c(0, 0))

  }

  if (legend.position.to.use[1] != "default") {
    graph.to.make <- graph.to.make +
      theme(legend.position = legend.position.to.use,
            legend.title = element_blank())
  }

  #### Add labels and title ####
  graph.to.make <- graph.to.make +
    labs(x = xlab.to.use,
         y = ylab.to.use,
         title = titleToUse)

  graph.to.make <- graph.to.make +
    theme(axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"))

  graph.to.make +
    coord_flip()

}
