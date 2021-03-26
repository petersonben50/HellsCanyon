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
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=''


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
cd ~/HellsCanyon/reports
rm -f */*_autoBinning*
mkdir /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binning/manualBinning/anvioDBs/original_summaries
cd /home/GLBRCORG/bpeterson26/HellsCanyon/code/
chmod +x executables/generate_large_bin_clusters_anvio.sh
condor_submit submission/generate_large_bin_clusters_anvio.sub


##########################
# Search bins for hgcA
##########################
screen -S HCC_binning
conda activate anvio6.2
cd /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binning/manualBinning/anvioDBs
mkdir original_summaries/hgcA_search
echo -e "assembly\tbin" >> original_summaries/hgcA_search/original_hgcA_bin_list.txt

cat /home/GLBRCORG/bpeterson26/HellsCanyon/metadata/lists/assembly_list.txt | while read assembly
do
  if [ ! -e original_summaries/uncurated_bins_from_$assembly ]; then
    echo "Summarizing binning from" $assembly
    anvi-summarize -c $assembly.db \
                    -p $assembly.merged/PROFILE.db \
                    -C CONCOCT \
                    -o original_summaries/uncurated_bins_from_$assembly

    ls original_summaries/uncurated_bins_from_$assembly/bin_by_bin | sed 's/\///' \
        > original_summaries/$assembly\_original_bin_list.txt
    cat original_summaries/$assembly\_original_bin_list.txt | while read bin
    do
      if [ -s original_summaries/uncurated_bins_from_$assembly/bin_by_bin/$bin/$bin-hgcaAnvio-hmm-sequences.txt ]; then
        echo $assembly$'\t'$bin >> original_summaries/hgcA_search/original_hgcA_bin_list.txt
      fi
    done
  else
    echo $assembly "already summarized"
  fi
done



####################################################
####################################################
# Run automatic binning algorithms
####################################################
####################################################

screen -S HCC_auto_binning
mkdir ~/HellsCanyon/dataEdited/binning/autoBinning
cd ~/HellsCanyon/dataEdited/binning/autoBinning
mkdir metabat2 maxbin2 dasTool finalBins
cd ~/HellsCanyon/reports
rm -f */*_autoBinning*

cd /home/GLBRCORG/bpeterson26/HellsCanyon/code/
chmod +x executables/automatic_binning.sh
condor_submit submission/automatic_binning.sub

cd /home/GLBRCORG/bpeterson26/HellsCanyon/metadata/lists
rm -f assembly_list.txt


####################################################
####################################################
# Add automatic binning information into anvio databases
####################################################
####################################################

screen -S HCC_anvioDBs_add_autoBins
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""

# Set up variables
metabat2=~/HellsCanyon/dataEdited/binning/autoBinning/metabat2
maxbin2=~/HellsCanyon/dataEdited/binning/autoBinning/maxbin2
dasTool=~/HellsCanyon/dataEdited/binning/autoBinning/finalBins
anvioDB=~/HellsCanyon/dataEdited/binning/manualBinning/anvioDBs
mkdir $anvioDB/S2B_files
mkdir $anvioDB/binning_collections

# Add in automatic binning data to anvio DB.
cat /home/GLBRCORG/bpeterson26/HellsCanyon/metadata/lists/assembly_list.txt | while read assembly
do

  echo "Adding scaffold to bin files for" $assembly

  # Copy all S2B files to one folder
  cp $metabat2/$assembly\_metabat_S2B.tsv \
      $maxbin2/$assembly\_maxbin_S2B.tsv \
      $dasTool/$assembly\_dasTool_S2B.tsv \
      $anvioDB/S2B_files

  anvi-import-collection $anvioDB/S2B_files/$assembly\_metabat_S2B.tsv \
                           -c $anvioDB/$assembly.db \
                           -p $anvioDB/$assembly.merged/PROFILE.db \
                           -C metabat2 \
                           --contigs-mode
  anvi-import-collection $anvioDB/S2B_files/$assembly\_maxbin_S2B.tsv \
                          -c $anvioDB/$assembly.db \
                          -p $anvioDB/$assembly.merged/PROFILE.db \
                          -C maxbin2 \
                          --contigs-mode
  anvi-import-collection $anvioDB/S2B_files/$assembly\_dasTool_S2B.tsv \
                          -c $anvioDB/$assembly.db \
                          -p $anvioDB/$assembly.merged/PROFILE.db \
                          -C dasTool \
                          --contigs-mode
  anvi-script-merge-collections -c $anvioDB/$assembly.db \
                                -i $anvioDB/S2B_files/$assembly*S2B.tsv \
                                -o $anvioDB/binning_collections/$assembly\_collections.tsv
done
