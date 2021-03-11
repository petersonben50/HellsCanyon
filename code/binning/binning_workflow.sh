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
mkdir binning/binning_automatic
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

cd /home/GLBRCORG/bpeterson26/HellsCanyon/metadata/lists
rm -f assembly_list.txt
