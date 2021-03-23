#!/bin/sh

##########################
# code/binning/binning_workflow.sh
# Benjamin D. Peterson

# This workflow will generate our set of
# hgcA+ genomes, as well as a set of
# uncurated bins.
##########################

####################################################
####################################################
# Prepare scaffolds and mapping files
####################################################
####################################################


##########################
# Filter out short scaffolds
##########################

screen -S HCC_binning
mkdir ~/HellsCanyon/dataEdited/binning
cd ~/HellsCanyon/dataEdited
mkdir binning/scaffolds
export PATH=/home/GLBRCORG/bpeterson26/miniconda3/bin:$PATH
source activate anvio6.2
PYTHONPATH=""
IFS=$'\n'


awk -F ',' '{ print $1 }' ~/HellsCanyon/metadata/lists/assembly_key.csv | while read assembly
do
  if [ ! -e binning/scaffolds/$assembly\_filtered_scaffolds.fna ]; then
    echo "Processing" $assembly "scaffolds for binning"
    anvi-script-reformat-fasta assemblies/scaffolds/$assembly\_assembly.fna \
                              -o binning/scaffolds/$assembly\_filtered_scaffolds.fna \
                              -l 2000
  else
    echo $assembly "scaffolds already processed for binning"
  fi
done


##########################
# Map reads to filtered scaffolds
##########################
cd /home/GLBRCORG/bpeterson26/HellsCanyon/reports/
rm -f outs/*_binningMapping.out \
      errs/*_binningMapping.err \
      logs/*_binningMapping.log

cd /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binning
mkdir mapping
mkdir mapping/indices

cd /home/GLBRCORG/bpeterson26/HellsCanyon/metadata/lists
awk -F ',' '{ print $1 }' assembly_key.csv > assembly_list.txt

cd /home/GLBRCORG/bpeterson26/HellsCanyon/code/
chmod +x executables/binning_mapping.sh
condor_submit submission/binning_mapping.sub


####################################################
####################################################
# Prepare anvi'o databases for manual binning
####################################################
####################################################
mkdir /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binning/manualBinning
mkdir /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binning/manualBinning/anvioDBs

##########################
# Generate contig databases
##########################
cd /home/GLBRCORG/bpeterson26/HellsCanyon/code/
chmod +x executables/anvio_DB_prep.sh
condor_submit submission/anvio_DB_prep.sub


##########################
# Generate read profiles
##########################
cd /home/GLBRCORG/bpeterson26/HellsCanyon/code/
chmod +x executables/anvio_profiling.sh
condor_submit submission/anvio_profiling.sub


##########################
# Run CONCOCT binning
##########################
cd /home/GLBRCORG/bpeterson26/HellsCanyon/code/
chmod +x executables/generate_large_bin_clusters_anvio.sh
condor_submit submission/generate_large_bin_clusters_anvio.sub







####################################################
####################################################
# Run automatic binning algorithms
####################################################
####################################################

screen -S HCC_auto_binning
mkdir ~/HellsCanyon/dataEdited/binning/autoBinning
cd ~/HellsCanyon/dataEdited/binning/autoBinning
mkdir metabat2 maxbin2 dasTool
cd ~/HellsCanyon/reports
rm -f */*_autoBinning*

cd /home/GLBRCORG/bpeterson26/HellsCanyon/code/
chmod +x executables/automatic_binning.sh
condor_submit submission/automatic_binning.sub

cd /home/GLBRCORG/bpeterson26/HellsCanyon/metadata/lists
rm -f assembly_list.txt
