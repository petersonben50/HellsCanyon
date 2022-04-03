#### code/hgcA_analysis/MeHg_hgcA_corr_main_figure.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon/")
library(lubridate)
library(tidyverse)
cb.translator <- readRDS("references/colorblind_friendly_colors.rds")


#### Read in data ####
Hg.data <- read.csv("dataEdited/waterChemistry/Hg_data.csv") %>%
  group_by(date, RM, depth, constituent) %>%
  summarise(concentration = sum(concentration))
hgcA.data <- read.csv("dataEdited/hgcA_analysis/depth/hgcA_coverage.csv")


#### Summarize hgcA data at each depth ####
hgcA.data <- hgcA.data %>%
  filter(rep == TRUE) %>%
  filter(!(RM %in% c(305, 314, 318)),
         !(RM == 286 & year(date) == 2019)) %>%
  group_by(date, RM, depth, redoxClassification) %>%
  summarise(hgcA_coverage = sum(coverage)) %>%
  ungroup()


#### Combine data ####
all.data <- full_join(Hg.data,
                      hgcA.data)
rm(Hg.data,
   hgcA.data)


#### Set up vectors for Hg and hgcA depth plots ####
color.vector <- c(cb.translator["black"],
                  cb.translator["bluishgreen"])
names(color.vector) <- c("hgcA_coverage", "FMHG")
labels.vector <- c("hgcA coverage\n(per 100X SCG)",
                   "Dissolved MeHg\n(ng/L)")
names(labels.vector) <- c("hgcA_coverage", "FMHG")
points.vector <- c(16, 17)
names(points.vector) <- c("hgcA_coverage", "FMHG")


#### Prepare data for plotting ####
MeHg.hgcA.data <- all.data %>%
  filter(year(date) %in% c(2017, 2018),
         RM %in% c(286, 300)) %>%
  spread(key = constituent,
         value = concentration) %>%
  filter(!is.na(hgcA_coverage)) %>%
  gather(key = constituent,
         value = abundance,
         -c(1:4)) %>%
  filter(constituent %in% names(color.vector)) %>%
  mutate(date.site = paste(year(date), "-RM", RM,
                           sep = ""))

#### Set up vectors for MeHg vs hgcA scatterplot ####
redox.color.vector <- c(cb.translator["bluishgreen"],
                        cb.translator["orange"],
                        cb.translator["blue"])
names(redox.color.vector) <- c("oxic", "suboxic", "euxinic")
redox.year.vector <- c(16, 17, 2)
names(redox.year.vector) <- c("2017", "2018", "2019")


#### Set up data for scatterplot ####
MeHg.hgcA.scatterplot.data <- all.data %>%
  spread(key = constituent,
         value = concentration) %>%
  filter(!is.na(hgcA_coverage)) %>%
  gather(key = constituent,
         value = abundance,
         -c(1:4)) %>%
  filter(constituent %in% names(color.vector)) %>%
  mutate(date.site = paste(year(date), "-RM", RM,
                           sep = "")) %>%
  ungroup() %>%
  spread(key = constituent,
         value = abundance) %>%
  mutate(redoxClassification = fct_relevel(redoxClassification,
                                           names(redox.color.vector)))



#### Set hgcA coverage at 20m at RM286 in 2018 approximately to DL of hgcA coverage ####
# We have a couple of lines that are essentially non-detects for hgcA.
# I think they actually are below the effective DL of our analysis. Due to the 
# log transformation, this is skewing the data left, so we're going to remove
# these samples.
MeHg.hgcA.scatterplot.data <- MeHg.hgcA.scatterplot.data %>%
  filter(hgcA_coverage > 0.001)


#### Linear regression of points ####
mehg.hgcA.model <- lm(log(FMHG, 10) ~ log(hgcA_coverage, 10),
                      data = MeHg.hgcA.scatterplot.data)
summary(mehg.hgcA.model)
# Get p-value
f <- summary(mehg.hgcA.model)$fstatistic
p.value <- pf(f[1],f[2],f[3],lower.tail=F) %>% round(4)

summary(mehg.hgcA.model)$coefficients[2, 1]


#### Generate scatterplot ####
hgcA.MeHg.scatterplot <- MeHg.hgcA.scatterplot.data %>%
  ggplot(aes(x = log(hgcA_coverage, 10),
             y = log(FMHG, 10),
             color = redoxClassification,
             shape = as.character(year(date)))) +
  geom_point(aes(color = redoxClassification)) +
  geom_abline(slope = coef(mehg.hgcA.model)[[2]],
              intercept = coef(mehg.hgcA.model)[[1]]) +
  geom_label(x = 0.5, y = -0.5,
             label = paste("Adjusted r2 = ", round(summary(mehg.hgcA.model)$adj.r.squared, 2), "\n",
                           "p = ", p.value,
                           sep = ""),
             color = "black") +
  scale_color_manual(values = redox.color.vector,
                     name = "Redox status") +
  scale_shape_manual(values = redox.year.vector,
                     name = "Year") +
  ylim(c(-1.5, 0.75)) +
  xlab("Log of hgcA abundance") +
  ylab("Log of dissolved MeHg (ng/L)") +
  theme_classic() +
  theme(legend.position = c(0.22, 0.75),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.key = element_rect(fill = "transparent", colour = "black"))
hgcA.MeHg.scatterplot



hgcA.MeHg.scatterplot.shading <- MeHg.hgcA.scatterplot.data %>%
  ggplot(aes(x = log(hgcA_coverage, 10),
             y = log(FMHG, 10))) +
  geom_smooth(method = lm ,
              color = "black",
              fill = "grey75",
              se = TRUE,
              level = 0.98) +
  geom_label(x = 0.5, y = -0.5,
             label = paste("Adjusted r2 = ", round(summary(mehg.hgcA.model)$adj.r.squared, 2), "\n",
                           "p = ", p.value,
                           sep = ""),
             color = "black") +
  scale_color_manual(values = redox.color.vector,
                     name = "Redox status") +
  scale_shape_manual(values = redox.year.vector,
                     name = "Year") +
  xlab("Log of hgcA abundance") +
  ylab("Log of dissolved MeHg (ng/L)") +
  ylim(c(-1.5, 0.75)) +
  theme_classic() +
  theme(legend.position = c(0.22, 0.75),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.key = element_rect(fill = "transparent", colour = "black"))
hgcA.MeHg.scatterplot.shading



#### Save out scatterplot ####
pdf("results/manuscript_figures/MeHg_hgcA_corr_main_figure.pdf",
    width = 5,
    height = 4.5)
hgcA.MeHg.scatterplot
dev.off()

pdf("results/manuscript_figures/MeHg_hgcA_corr_main_figure_shading.pdf",
    width = 5,
    height = 4.5)
hgcA.MeHg.scatterplot.shading
dev.off()
