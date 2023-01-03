#### code/hgcA_analysis/MeHg_hgcA_corr_main_figure.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(tidyverse)
source("code/HCC_plotting_needs.R")


#### Read in data ####
geochem.data <- readRDS("dataEdited/geochem/geochem_WC_adjusted_for_redox_classification.rds")
hgcA.data <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv")


#### Summarize hgcA data at each depth ####
hgcA.data.sum <- hgcA.data %>%
  filter(rep == TRUE) %>%
  filter(!(RM %in% c(305, 314, 318)),
         !(RM == 286 & year(date) == 2019)) %>%
  group_by(date, RM, depth, redoxClassification) %>%
  summarise(hgcA_coverage = sum(coverage)) %>%
  mutate(RM = as.character(RM),
         date = as.Date(date)) %>%
  ungroup()


#### Combine data ####
all.data <- right_join(geochem.data,
                      hgcA.data.sum)
rm(geochem.data,
   hgcA.data.sum)



#### Set hgcA coverage at 20m at RM286 in 2018 approximately to DL of hgcA coverage ####
# We have a couple of lines that are essentially non-detects for hgcA.
# I think they actually are below the effective DL of our analysis. Due to the 
# log transformation, this is skewing the data left, so we're going to remove
# these samples.
all.data$hgcA_coverage[all.data$hgcA_coverage <= 0.001] <- 0.001



#### Quick look ####
all.data %>%
  ggplot(aes(x = log(hgcA_coverage, 10),
             y = log(MeHg_diss_ngL, 10))) +
  geom_smooth(method = lm ,
              color = "black",
              fill = "grey75",
              se = TRUE,
              fullrange = TRUE,
              level = 0.95) +
  geom_point(aes(color = redox_status,
                 shape = as.character(year(date))),
             size = 2)


#### Linear regression of points ####
mehg.hgcA.model <- lm(log(all.data$MeHg_diss_ngL, 10) ~ log(all.data$hgcA_coverage, 10),
                      data = all.data)
summary(mehg.hgcA.model)
# Get p-value
f <- summary(mehg.hgcA.model)$fstatistic
p.value <- pf(f[1],f[2],f[3],lower.tail=F) %>% round(4)

summary(mehg.hgcA.model)$coefficients[2, 1]


#### Generate scatterplot ####
hgcA.MeHg.scatterplot <- all.data %>%
  ggplot(aes(x = hgcA_coverage,
             y = MeHg_diss_ngL)) +
  geom_smooth(method = lm ,
              color = "black",
              fill = "grey75",
              se = TRUE,
              fullrange = TRUE,
              level = 0.95) +
  geom_point(aes(color = redox_status,
                 shape = as.character(year(date))),
             size = 2) +
  geom_abline(slope = coef(mehg.hgcA.model)[[2]],
              intercept = coef(mehg.hgcA.model)[[1]]) +
  geom_label(x = 0.5, y = -0.5,
             label = paste("Adjusted R2 = ", round(summary(mehg.hgcA.model)$adj.r.squared, 2), "\n",
                           "p = ", p.value,
                           sep = ""),
             color = "black") +
  scale_color_manual(values = color.vector,
                     labels = renaming.vector,
                     name = "Redox status") +
  scale_shape_manual(values = shape.vector,
                     labels = renaming.vector,
                     name = "Year") +
  xlab("hgcA abundance (%)") +
  ylab("Filter-passing MeHg (ng/L)") +
  theme_classic() +
  theme(legend.position = c(0.22, 0.75),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.key = element_rect(fill = "transparent", colour = "black")) +
  scale_x_continuous(limits = c(0.001, 10),
                     trans = 'log10') +
  scale_y_continuous(limits = c(0.05, 4),
                     trans = 'log10')
hgcA.MeHg.scatterplot



#### Save out scatterplot ####
pdf("results/hgcA_analysis/MeHg_vs_hgcA.pdf",
    width = 3.5,
    height = 3)
hgcA.MeHg.scatterplot
dev.off()
