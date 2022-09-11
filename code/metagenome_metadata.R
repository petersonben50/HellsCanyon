#### code/metagenome_metadata.R ####
# Benjamin D. Peterson

# This file contains the code to generate nice tables with
# information on the metagenome sequencing and assembly.

#### Always start with a clean slate, ya bum ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/HellsCanyon")
library(dplyr)
library(gridExtra)
library(readxl)


#### Generate mapping key ####
MG.list <- read.csv("metadata/lists/metagenome_key.csv",
                    header = FALSE,
                    col.names = c("metagenomeID", "samplingYear"))
assembly.list <- read.csv("metadata/lists/assembly_key.csv",
                          header = FALSE,
                          col.names = c("assemblyID", "samplingYear"))
mapping.key <- full_join(MG.list,
                         assembly.list) %>%
  select(metagenomeID, assemblyID, samplingYear)

write.table(mapping.key,
            "metadata/lists/mapping_key.tsv",
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)


#### Generate table of MG sample site information ####

# Read in metagenome metadata
MG.metadata <- read_xlsx("dataEdited/dnaSequencing/metagenome_key.xlsx")

# Read in extraction metadata
extraction.metadata <- read_xlsx("dataEdited/dnaExtractions/DNA_extractions_data.xlsx")

# Read in filter metadata
filter.metadata <- read_xlsx("metadata/NA_filters.xlsx")

all.metadata <- left_join(MG.metadata,
                          extraction.metadata,
                          by = "extractionID") %>%
  left_join(filter.metadata,
            by = "filterID")



#### Save out metadata for analysis ####

all.metadata %>%
  select(metagenomeID, date, RM, depth) %>%
  mutate(RM = gsub("RM", "", RM))%>%
  write.csv("metadata/metagenome_metadata.csv",
            row.names = FALSE)


#### BioSample entry prep ####
lat.lon.vector <- c("44.8229 N 116.9156 W",
                    "44.6934 N 117.0778 W",
                    "44.6289 N 117.1066 W",
                    "44.5710 N 117.1443 W",
                    "44.5244 N 117.1717 W",
                    "44.48708 N 117.2171 W")
names(lat.lon.vector) <- c("RM286", "RM300", "RM305", "RM310", "RM314", "RM318")
sampleIDs.of.interest <- unique(all.metadata$filterID)
all.metadata %>%
  mutate(
    # Use descriptive name for sample name
    sample_name = paste("HellsCanyonBrownleeRes", "_",
                        year(date), "_",
                        month(date, label = TRUE, abbr = TRUE), "_",
                        RM, "_",
                        as.character(depth), "m",
                        sep = ""),
    # Set sample title as sample ID from project
    sample_title = filterID,
    bioproject_accession = "PRJNA878929",
    # "Organism" name selected from this list: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Tree&id=410657&lvl=3&lin=f&keep=1&srchmode=1&unlock
    # Instructions to do so here: https://www.ncbi.nlm.nih.gov/biosample/docs/organism/#metagenomes
    organism = "freshwater metagenome",
    collection_date = date,
    # env_broad_scale selected from this list: https://ontobee.org/ontology/ENVO?iri=http://purl.obolibrary.org/obo/ENVO_00000428
    # Instructions to do so here: https://www.gensc.org/pages/standards/all-terms.html
    env_broad_scale = "lake [ENVO:00000020]",
    # Same as above
    env_local_scale = "hydroelectric reservoir [ENVO:00000448]|lake with an anoxic hypolimnion [ENVO:01001073]|freshwater lake [ENVO:00000021]",
    # Same as above
    env_medium = "fresh water [ENVO_00002011]",
    geo_loc_name = "USA: Hells Canyon Complex - Brownlee Reservoir",
    lat_lon = lat.lon.vector[RM]) %>%
  select(sample_name, sample_title, bioproject_accession, organism, collection_date,
         depth, env_broad_scale, env_local_scale, env_medium, geo_loc_name, lat_lon) %>%
  write.csv("dataEdited/metagenomes/NCBI_upload/temp_NCBI_info_MIMS.csv",
            quote = FALSE,
            row.names = FALSE)



#### SRA entry prep ####
biosample.ids <- read.table("dataEdited/metagenomes/NCBI_upload/NCBI_info_MIMS_with_BioSample.txt",
                            sep = '\t',
                            skip = 1) %>%
  select(V1, V3)
names(biosample.ids) <- c("biosample_accession", "sample_name")
vector.for.prep <- c(rep("Library prep and sequencing done at QB3 Genomics Center at Berkeley. Library prep used Kapa Biosystem Library Prep kit with a target insert length of ∼600 bp. Sequenced using the S4 chemistry. 150 bp paired-end sequencing", 20),
                     rep("Library prep and sequencing done at Michigan State University Genomics Core. Library prep used TakaraBio SMARTer ThruPLEX DNA. 150 bp paired-end sequencing", 14))
names(vector.for.prep) <- all.metadata$metagenomeID
vector.for.instrument <- c(rep("Illumina NovaSeq 6000", 20),
                           rep("Illumina HiSeq 4000", 14))
names(vector.for.instrument) <- all.metadata$metagenomeID


all.metadata %>%
  mutate(
    # Use descriptive name for sample name
    sample_name = paste("HellsCanyonBrownleeRes", "_",
                        year(date), "_",
                        month(date, label = TRUE, abbr = TRUE), "_",
                        RM, "_",
                        as.character(depth), "m",
                        sep = "")) %>%
    left_join(biosample.ids) %>%
  rename(library_ID = metagenomeID) %>%
  mutate(title = paste("Metagenomic sequencing of microbial community from the water column of Brownlee Reservoir on ",
                       date, " at ", RM, " from ", depth, "m",
                       sep = ""),
         library_strategy = "WGS",
         library_source = "METAGENOMIC",
         library_selection = "RANDOM",
         library_layout = "paired",
         platform = "ILLUMINA",
         instrument_model = vector.for.instrument[library_ID],
         design_description = vector.for.prep[library_ID],
         filetype = "fastq",
         filename = paste(library_ID, "_R1.fastq.gz",
                          sep = ""),
         filename2 = paste(library_ID, "_R2.fastq.gz",
                           sep = ""),
         filename3 = "",
         filename4 = "",
         assembly = "",
         fasta_file = "") %>%
  select(biosample_accession, library_ID, title, library_strategy, library_source,
         library_selection, library_layout, platform, instrument_model, design_description,
         filetype, filename, filename2, filename3, filename4, assembly, fasta_file) %>%
  arrange(biosample_accession) %>%
  write.csv(file = "dataEdited/metagenomes/NCBI_upload/temp_NCBI_info_SRA.csv",
            quote = FALSE,
            row.names = FALSE)


