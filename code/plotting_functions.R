#### code/gene_plotting_functions.R ####
# Written by Benjamin D. Peterson



#### Set redox plotting function for TEAP figure ####
plot.DO.nitrate.Mn.sulfide.profile <- function(geochem.data.of.interest = geochem.data,
                                               RM.of.interest,
                                               date.of.interest,
                                               depth.range.of.interest,
                                               concentration.of.interest,
                                               legend.location.of.interest,
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

  geochem.data.of.interest %>%
    ggplot(aes(x = depth,
               y = concentration,
               color = constituent,
               linetype = constituent,
               shape = constituent)) +
    geom_point(aes(color = constituent),
               size = 2) +
    geom_line(aes(color = constituent)) +
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
}



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
                                            legend.position.to.use = "default",
                                            legend.title.to.use = element_blank(),
                                            titleToUse = NULL) {

  marker.depth.df[, "gene.names.to.use.for.filtering"] <- marker.depth.df[, gene.name.column]

  clean.coverage <- marker.depth.df %>%
    filter(gene.names.to.use.for.filtering %in% genesOfInterest) %>%
    filter(year(ymd(date)) == yearOfInterest) %>%
    filter(RM == RMofInterest) %>%
    group_by(gene.names.to.use.for.filtering, depth) %>%
    summarise(coverage = sum(coverage)) %>%
    ungroup() %>%
    arrange(depth)

  graph.of.interest <- clean.coverage %>%
    ggplot(aes(x = depth,
               y = coverage)) +
    theme_classic() +
    scale_x_reverse(limits = depth_limits)  +
    theme(axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))

  #### Add labels and title ####
  if (is.null(titleToUse)) {
    titleToUse = paste("Coverage of",
                       paste(genesOfInterest, collapse = ","))

  }
  graph.of.interest <- graph.of.interest +
    labs(title = titleToUse,
         y = xlab.to.use,
         x = "Depth (m)")

  #### Add colors if defined in the call ####
  if (!is.null(color.vector.to.use)) {
    graph.of.interest <- graph.of.interest +
      geom_line(aes(group = gene.names.to.use.for.filtering,
                    color = gene.names.to.use.for.filtering,
                    linetype = gene.names.to.use.for.filtering)) +
      geom_point(aes(group = gene.names.to.use.for.filtering,
                     color = gene.names.to.use.for.filtering,
                     shape = gene.names.to.use.for.filtering)) +
      scale_colour_manual(values=color.vector.to.use)
  } else {
    graph.of.interest <- graph.of.interest +
      geom_point(shape = gene.names.to.use.for.filtering) +
      geom_line(aes(group = gene.names.to.use.for.filtering,
                    linetype = gene.names.to.use.for.filtering))
  }

  #### Add point shapes if defined in the call ####
  if (!is.null(point.vector.to.use)) {
    graph.of.interest <- graph.of.interest +
      scale_shape_manual(values = point.vector.to.use)
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
      scale_y_continuous(limits = c(0, max.coverage))

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
            legend.title = legend.title.to.use,
            legend.box.background = element_rect(colour = "black"))
  }

  graph.of.interest + coord_flip()

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
