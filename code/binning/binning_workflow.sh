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


################################################
################################################
# Manually bin hgcA+ bins
################################################
################################################

screen -S HCC_anvioDBs_binning
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""
cd ~/HellsCanyon/dataEdited/binning/manualBinning/
# Copy the database folder
#cp -avr anvioDBs anvioDBs_modified
#cat /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binning/manualBinning/anvioDBs/original_summaries/hgcA_search/original_hgcA_bin_list.txt
assembly=fall2017coassembly
bin=Bin_77
anvi-refine -p anvioDBs_modified/$assembly.merged/PROFILE.db \
            -c anvioDBs_modified/$assembly.db \
            -C CONCOCT \
            -b $bin \
            -A anvioDBs_modified/binning_collections/$assembly\_collections.tsv \
            --taxonomic-level "t_phylum"




################################################
################################################
# Summarize/export curated bins
################################################
################################################

##########################
# Rename and summarize bins
##########################
screen -S HCC_binsPostProcessing
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""
PERL5LIB=""
cd ~/HellsCanyon/dataEdited/binning/manualBinning/anvioDBs_modified

# Summarize them
cat /home/GLBRCORG/bpeterson26/HellsCanyon/metadata/lists/assembly_list.txt | while read assembly
do
  # If old summary exists that we want to delete, uncomment
  # the following:
  # rm -r $assembly.summary.curated
  if [ ! -d $assembly.curated.summary ]; then
    echo "Summarizing bins for" $assembly
    anvi-summarize -c $assembly.db \
                    -p $assembly.merged/PROFILE.db \
                    -C CONCOCT \
                    -o $assembly.curated.summary
  else
    echo "We already summarized the curated bins for" $assembly
  fi
done
conda deactivate


##########################
# Pull out DNA files from hgcA+ bins
##########################
# First need to set up new directory
cd ~/HellsCanyon/dataEdited/binning/manualBinning/
mkdir binsRaw
mkdir binsRaw/DNA
binsRaw=~/HellsCanyon/dataEdited/binning/manualBinning/binsRaw

cat /home/GLBRCORG/bpeterson26/HellsCanyon/metadata/lists/assembly_list.txt | while read assembly
do
  binSummary=~/HellsCanyon/dataEdited/binning/manualBinning/anvioDBs_modified/$assembly.curated.summary
  if [ -e $binSummary ]; then
    if [ ! -e $binsRaw/DNA/$assembly* ]; then
      cd $binSummary/bin_by_bin
      ls | sed 's/\///' | while read bin
      do
        isThereHgcA=`cat $bin/$bin\-hgcaAnvio-hmm-sequences.txt | wc -l`
        if [ ! $isThereHgcA -eq 0 ]; then
          echo "Copying" $bin "to binsRaw folder"
          cp $bin/$bin-contigs.fa $binsRaw/DNA/$bin.fna
        else
          echo "No hgcA in" $bin
        fi
      done
    else
      echo "Hey, there are some bins from" $assembly "already in here"
      echo "You might wanna check that out before you start overwriting stuff"
    fi
  else
    echo "Summarize anvioDB for" $assembly", dummy."
  fi
done

# Generate list of hgcA+ bins
cd $binsRaw/DNA
ls *.fna | \
  sed 's/.fna//' \
  > binsRaw_hgcA_list.txt






####################################################
####################################################
# Check quality of bins
####################################################
####################################################

##########################
# Completeness/redundancy estimates from anvio
##########################
binsRaw=~/HellsCanyon/dataEdited/binning/manualBinning/binsRaw
mkdir $binsRaw/anvio_data
mkdir $binsRaw/anvio_data/completeness_redundancy

# Copy summary files into a single folder.
cat /home/GLBRCORG/bpeterson26/HellsCanyon/metadata/lists/assembly_list.txt | while read assembly
do
  binSummary=~/HellsCanyon/dataEdited/binning/manualBinning/anvioDBs_modified/$assembly.curated.summary
  if [ -e $binSummary/bins_summary.txt ]; then
    cp $binSummary/bins_summary.txt $binsRaw/anvio_data/completeness_redundancy/$assembly\_bins_summary.txt
  else
    echo $assembly "has not been summarized."
  fi
done

# Concatenate summaries into a single file.
cd $binsRaw/anvio_data/completeness_redundancy
head -n 1 fall2017cluster1_bins_summary.txt > bins_summary_all.txt
ls *_bins_summary.txt | while read file
do
  tail -n +2 $file >> bins_summary_all.txt
done

# Only keep the summaries for the hgcA+ bins.
head -n 1 bins_summary_all.txt > bins_summary_hgcA.txt
cat $binsRaw/DNA/binsRaw_hgcA_list.txt | while read hgcA_bin
do
  grep $hgcA_bin bins_summary_all.txt >> bins_summary_hgcA.txt
done

head -n 1 bins_summary_hgcA.txt > bins_summary_hgcA_good.txt
awk '{ if (($7 > 50) && ($8 < 10)) print $0 }' bins_summary_hgcA.txt >> bins_summary_hgcA_good.txt
tail -n +2 bins_summary_hgcA_good.txt | \
  awk '{ print $1 }' \
  > bins_list_hgcA_good.txt




##########################
# Completeness/redundancy estimates from CheckM
##########################

screen -S HCC_checkM
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
conda activate bioinformatics
binsRaw=~/HellsCanyon/dataEdited/binning/manualBinning/binsRaw

cd $binsRaw
if [ -d checkM ]; then
  echo "Removing old checkM folder"
  rm -rf checkM
fi
mkdir checkM
checkm lineage_wf \
      -x .fna \
      -t 16 \
      DNA \
      checkM
checkm qa checkM/lineage.ms \
        checkM \
        -o 2 \
        -f checkM/checkm.out \
        --tab_table
awk -F '\t' \
  -v OFS=',' \
  '{ print $1,$6,$7,$8,$9,$11,$13,$15,$17,$19,$23 }' \
  checkM/checkm.out \
  > checkM/checkM_stats.csv

# Download checkM/checkM_stats.csv to local computer:
# dataEdited/2017_analysis_bins/binning/rawBins/bin_quality
cd ~/HellsCanyon/dataEdited/binning/manualBinning
mkdir binsGood
mkdir binsGood/DNA
mkdir binsGood/checkM
awk -F ',' '{ if (($2 > 50) && ($3 < 10)) print $0 }' \
  binsRaw/checkM/checkM_stats.csv \
  > binsGood/checkM/good_bins_data.txt
awk -F ',' '{ print $1 }' binsGood/checkM/good_bins_data.txt \
  > binsGood/checkM/good_bins_list.txt
cat binsGood/checkM/good_bins_list.txt | while read binsGood
do
  echo "Copying" $binsGood
  cp binsRaw/DNA/$binsGood.fna binsGood/DNA
done

cd binsGood/DNA
scripts=~/Everglades/code/generalUse
ls *fna | while read fna
do
  echo "Cleaning" $fna
  python $scripts/cleanFASTA.py $fna
  mv -f $fna\_temp.fasta $fna
done


####################################################
####################################################
# Check out taxonomy of bins with GTDB
####################################################
####################################################

screen -S HCC_GTDB
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
conda activate gtdbtk
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood
rm -rf taxonomy
mkdir taxonomy

gtdbtk classify_wf \
        --cpus 16 \
        --extension fna \
        --genome_dir ./DNA \
        --out_dir taxonomy
# Summarize them
cd taxonomy
grep -h 'anvio_hgcA' gtdbtk.*.summary.tsv \
        | awk -F '\t' '{ print $1"\t"$2 }' \
        > taxonomy_summary.txt
conda deactivate
