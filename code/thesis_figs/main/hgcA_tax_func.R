#### code/hgcA_analysis/hgcA_abundance.R ####
# Written by Benjamin D. Peterson


#### Get ready. Get set ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(lubridate)
library(patchwork)
library(readxl)
library(tidyverse)
source("code/gene_plotting_functions.R")
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
cb.translator <- c(cb.translator, "gray80")
names(cb.translator)[length(cb.translator)] <- "gray80"
cb.translator <- c(cb.translator, "gray50")
names(cb.translator)[length(cb.translator)] <- "gray50"
cb.translator <- c(cb.translator, "gray20")
names(cb.translator)[length(cb.translator)] <- "gray20"

#### Read in hgcA classification ####
tax.data <- read_xlsx("dataEdited/hgcA_analysis/hgcA_information_edited.xlsx") %>%
  filter(usedForAbundance == TRUE) %>%
  select(seqID, manual_classification, predicted_metabolism)



#### Read in depth data ####
hgcA.data <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv") %>%
  filter(rep) %>%
  left_join(tax.data)


#### Make color vector for taxonomy ####
hgcA.manual.taxonomy <- read_xlsx("dataEdited/hgcA_analysis/phylogeny/manual_taxonomy.xlsx",
                                  sheet = "colors_to_use")
color.vector.tax <- cb.translator[hgcA.manual.taxonomy$colorsToUse]
names(color.vector.tax) <- hgcA.manual.taxonomy$seqID


#### Make color vector for function ####
hgcA.manual.function <- read_xlsx("dataEdited/hgcA_analysis/hgcA_information_edited.xlsx",
                                  sheet = "metabolism_colors_to_use")
color.vector.fun <- cb.translator[hgcA.manual.function$colorsToUse]
names(color.vector.fun) <- hgcA.manual.function$metabolismID


#### Combine rare taxa into "other" ####
hgcA.data$manual_classification[which(hgcA.data$manual_classification %in% c("Actinobacteria",
                                                                             "Spirochaetes"))] <- "other"



#### Generate plots for 2017 at RM286 ####
hgcA.286.2017 <- plot.profile.of.gene.with.taxonomy.function.noDepth(data.to.use = hgcA.data,
                                                                     RM.of.interest = 286,
                                                                     year.of.interest = 2017,
                                                                     legend.position.to.use = "none",
                                                                     ylab.to.use = "Gene abundance)")
pdf("results/manuscript_figures/hgcA_tax_fig/286_2017.pdf",
    width = 2.2,
    height = 3.25)
hgcA.286.2017
dev.off()




#### Generate zoomed-in plots for 2017 at RM286 ####
hgcA.data.truncated <- hgcA.data %>%
  filter(depth < 60)
hgcA.286.2017.zoomIn <- plot.profile.of.gene.with.taxonomy.function.noDepth(data.to.use = hgcA.data.truncated,
                                                                            RM.of.interest = 286,
                                                                            year.of.interest = 2017,
                                                                            legend.position.to.use = "none",
                                                                            ylab.to.use = element_blank(),
                                                                            xlab.to.use = element_blank())
pdf("results/manuscript_figures/hgcA_tax_fig/286_2017_zoom.pdf",
    width = 1.3,
    height = 2.2)
hgcA.286.2017.zoomIn
dev.off()




#### Generate plots for 2017 at RM300 ####
hgcA.300.2017 <- plot.profile.of.gene.with.taxonomy.function.noDepth(data.to.use = hgcA.data,
                                                                     RM.of.interest = 300,
                                                                     year.of.interest = 2017,
                                                                     legend.position.to.use = "none",
                                                                     ylab.to.use = "Gene abundance")
pdf("results/manuscript_figures/hgcA_tax_fig/300_2017.pdf",
    width = 2,
    height = 2)
hgcA.300.2017
dev.off()



#### Generate plots for 2018 at RM286 ####
hgcA.286.2018 <- plot.profile.of.gene.with.taxonomy.function.noDepth(data.to.use = hgcA.data,
                                                                     RM.of.interest = 286,
                                                                     year.of.interest = 2018,
                                                                     legend.position.to.use = "none",
                                                                     ylab.to.use = "Gene abundance")
pdf("results/manuscript_figures/hgcA_tax_fig/286_2018.pdf",
    width = 2.2,
    height = 3.25)
hgcA.286.2018
dev.off()


#### Generate zoomed-in plots on upper water column for 2018 at RM286 ####
hgcA.data.lower <- hgcA.data %>%
  filter(depth >= 38)
hgcA.286.2018.lower <- plot.profile.of.gene.with.taxonomy.function.noDepth(data.to.use = hgcA.data.lower,
                                                                           RM.of.interest = 286,
                                                                           year.of.interest = 2018,
                                                                           legend.position.to.use = "none",
                                                                           ylab.to.use = element_blank(),
                                                                           xlab.to.use = element_blank())
pdf("results/manuscript_figures/hgcA_tax_fig/286_2018_lower.pdf",
    width = 2,
    height = 2)
hgcA.286.2018.lower
dev.off()



#### Generate plots for 2018 at RM300 ####
hgcA.300.2018 <- plot.profile.of.gene.with.taxonomy.function.noDepth(data.to.use = hgcA.data,
                                                                     RM.of.interest = 300,
                                                                     year.of.interest = 2018,
                                                                     legend.position.to.use = "none",
                                                                     ylab.to.use = "Gene abundance")
pdf("results/manuscript_figures/hgcA_tax_fig/300_2018.pdf",
    width = 2,
    height = 2)
hgcA.300.2018
dev.off()


#### Generate plots for 2019 at RM300 ####
hgcA.300.2019 <- plot.profile.of.gene.with.taxonomy.function.noDepth(data.to.use = hgcA.data,
                                                                     RM.of.interest = 300,
                                                                     year.of.interest = 2019,
                                                                     legend.position.to.use = "none",
                                                                     ylab.to.use = "Gene abundance")
pdf("results/manuscript_figures/hgcA_tax_fig/300_2019.pdf",
    width = 2,
    height = 1.33)
hgcA.300.2019
dev.off()




#### Generate plots for 2019 at RM310 ####
hgcA.310.2019 <- plot.profile.of.gene.with.taxonomy.function.noDepth(data.to.use = hgcA.data,
                                                                     RM.of.interest = 310,
                                                                     year.of.interest = 2019,
                                                                     legend.position.to.use = "none",
                                                                     ylab.to.use = "Gene abundance")
pdf("results/manuscript_figures/hgcA_tax_fig/310_2019.pdf",
    width = 2,
    height = 1.33)
hgcA.310.2019
dev.off()



#### Get plot of legend ####

tax.order.by.function <- hgcA.data %>%
  arrange(predicted_metabolism) %>%
  select(manual_classification) %>%
  unlist(use.names = FALSE) %>%
  unique() %>%
  as.character()

hgcA.data <- hgcA.data %>%
  mutate(manual_classification = fct_relevel(manual_classification,
                                             tax.order.by.function))

plot.for.legend <- ggplot(hgcA.data,
                          aes(x = depth,
                              y = coverage,
                              fill = manual_classification)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color.vector.tax)


pdf("results/manuscript_figures/hgcA_tax_fig/legend_to_use.pdf",
    width = 2,
    height = 2.75)
grid.newpage()
grid.draw(cowplot::get_legend(plot.for.legend))
dev.off()

