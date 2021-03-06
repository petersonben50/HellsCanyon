#### code/SCGs_depth_cleaning.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/HellsCanyon")
library(tidyverse)


#### Read in data ####
scg.abundance <- read.table("dataEdited/scg_abundance/scg_coverage.tsv",
                            sep = "\t",
                            stringsAsFactors = FALSE,
                            col.names = c("geneID", "metagenomeID", "coverage")) 



#### Plot out raw data ####
pdf("results/scg_abundance/scg_abundance.pdf",
    height = 4.5,
    width = 6)
scg.abundance %>%
  ggplot(aes(x = metagenomeID,
             y = coverage)) +
  geom_line(mapping = aes(group = geneID)) +
  geom_point() +
  theme_classic() +
  stat_summary(geom = "point", fun = "mean",
               col = "black", fill = "red",
               size = 3, shape = 24) +
  theme(axis.text.x = element_text(angle = 90))
dev.off()


#### Generate normalization factor for each metagenomes ####
mean.scg.abundance <- scg.abundance %>%
  group_by(metagenomeID) %>%
  summarize(coverage = mean(coverage))
# Normalize to a coverage of 100
normalized.mean.scg.abundance <- mean.scg.abundance %>%
  mutate(NF = 100 / coverage)


#### Generate normalization vector for each metagenomes ####

normalization.vector <- normalized.mean.scg.abundance$NF
names(normalization.vector) <- normalized.mean.scg.abundance$metagenomeID
saveRDS(normalization.vector,
        "dataEdited/scg_abundance/scg_normalization_vector.rds")
write.csv(normalized.mean.scg.abundance,
          "results/scg_abundance/scg_abundance_normalized_table.csv",
          row.names = FALSE)
