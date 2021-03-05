#!/bin/sh

######################
# code/SCG_coverage_for_normalization.sh
# Benjamin D. Peterson
######################

######################
# Get set up
######################
cd ~/HellsCanyon/metadata/lists
awk -F ',' '{ print $1 }' assembly_key.csv > assembly_list.txt
awk -F ',' '{ print $1 }' metagenome_key.csv > metagenome_list.txt
mkdir /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/scg_abundance

######################
# Submit jobs
######################
chmod +x /home/GLBRCORG/bpeterson26/HellsCanyon/code/executables/SCG_abundance_in_assemblies.sh
condor_submit /home/GLBRCORG/bpeterson26/HellsCanyon/code/submission/SCG_abundance_in_assemblies.sub

######################
# Clean up
######################
cd ~/HellsCanyon/metadata/lists
rm -f assembly_list.txt metagenome_list.txt
