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

# Move all bins from Das Tool into final bin set
cd ~/HellsCanyon/dataEdited/binning/autoBinning
mv finalBins dasToolOutput
mkdir completeBinSet
mkdir completeBinSet/DNA
assembly_list=/home/GLBRCORG/bpeterson26/HellsCanyon/metadata/lists/assembly_list.txt
cat $assembly_list | while read assemblyID
do
  echo "Moving bins from" $assemblyID
  cp dasToolOutput/$assemblyID\_bins/*fna completeBinSet/DNA/
done



####################################################
####################################################
# Initial bin characterization
####################################################
####################################################

##########################
# Get quality bins
##########################

screen -S HCC_checkM
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
conda activate bioinformatics
binsRaw=~/HellsCanyon/dataEdited/binning/autoBinning

cd $binsRaw
if [ -d checkM ]; then
  echo "Removing old checkM folder"
  rm -rf checkM
fi
mkdir checkM
checkm lineage_wf \
      -x .fna \
      -t 16 \
      completeBinSet/DNA \
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
cd ~/HellsCanyon/dataEdited/binning/autoBinning
mkdir hqBinSet
mkdir hqBinSet/DNA
mkdir hqBinSet/checkM
awk -F ',' '{ if (($2 > 50) && ($3 < 10)) print $0 }' \
  checkM/checkM_stats.csv \
  > hqBinSet/checkM/good_bins_data.txt
awk -F ',' '{ print $1 }' hqBinSet/checkM/good_bins_data.txt \
  > hqBinSet/checkM/good_bins_list.txt

scripts=~/HellsCanyon/code/generalUse
cat hqBinSet/checkM/good_bins_list.txt | while read binID
do
  if [ -e completeBinSet/DNA/$binID.fna ]; then
    if [ -e hqBinSet/DNA/$binID.fna ]; then
      echo "Already copied" $binID
    else
      echo "Copying over" $binID.fna
      cp completeBinSet/DNA/$binID.fna hqBinSet/DNA/$binID.fna
      python $scripts/cleanFASTA.py hqBinSet/DNA/$binID.fna
      mv -f hqBinSet/DNA/$binID.fna_temp.fasta hqBinSet/DNA/$binID.fna
    fi
  else
    echo $binID "doesn't exist"
  fi
done


##########################
# Run GTDB
##########################

screen -S HCC_GTDB
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
conda activate gtdbtk
cd ~/HellsCanyon/dataEdited/binning/autoBinning/hqBinSet
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




##########################
# Run ORF prediction on all bins
##########################
reports=/home/GLBRCORG/bpeterson26/HellsCanyon/reports
mkdir $reports/outs/orfPredictionBins $reports/errs/orfPredictionBins $reports/logs/orfPredictionBins
cd /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binning/autoBinning/hqBinSet
mkdir ORFs

cd /home/GLBRCORG/bpeterson26/HellsCanyon/code/
chmod +x executables/ORF_prediction_bins.sh
condor_submit submission/ORF_prediction_bins.sub



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


##########################
# Final bin list
##########################
cd ~/HellsCanyon/dataEdited/binning/manualBinning
cat binsGood/checkM/good_bins_list.txt binsRaw/anvio_data/completeness_redundancy/bins_list_hgcA_good.txt | \
  sort | uniq \
  > goodBins_list.txt


cat goodBins_list.txt | while read binsGood
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



####################################################
####################################################
# Pull out coverage information from anvio
####################################################
####################################################

# Coverage data from anvio
binsGood=~/HellsCanyon/dataEdited/binning/manualBinning/binsGood
mkdir $binsGood/coverageAnvio
cd $binsGood/coverageAnvio

cat /home/GLBRCORG/bpeterson26/HellsCanyon/metadata/lists/assembly_list.txt | while read assembly
do
  summary=~/HellsCanyon/dataEdited/binning/manualBinning/anvioDBs_modified/$assembly.curated.summary/bins_across_samples
  if [ -e $summary/mean_coverage_Q2Q3.txt ]; then
    if [ ! -e $binsGood/coverageAnvio/$assembly\_coverage.txt ]; then
      cp -f $summary/mean_coverage_Q2Q3.txt $binsGood/coverageAnvio/$assembly\_coverage.txt
    fi
  else
    echo $assembly "has not been summarized."
  fi
done


awk -F ',' '{ print $2 }' ~/HellsCanyon/metadata/lists/assembly_key.csv | \
  sort | uniq | \
  while read year
do
  initialAssembly=`awk -F ',' -v year="$year" '$2 == year { print $1 }' ~/HellsCanyon/metadata/lists/assembly_key.csv | head -n 1`
  head -n 1 $initialAssembly\_coverage.txt > coverage_$year.txt

  awk -F ',' -v year="$year" '$2 == year { print $1 }' ~/HellsCanyon/metadata/lists/assembly_key.csv | while read assembly
  do
    tail -n +2 $assembly\_coverage.txt >> coverage_$year.txt
  done

  head -n 1 coverage_$year.txt > coverage_goodBins_$year.txt
  grep -f ~/HellsCanyon/dataEdited/binning/manualBinning/goodBins_list.txt coverage_$year.txt >> coverage_goodBins_$year.txt
  rm -f coverage_$year.txt
done
rm -f *_coverage.txt
# Download the "coverage_goodBins_201*.txt" files to my computer



####################################################
####################################################
# Get ORFs for bins
####################################################
####################################################

screen -S HCC_binsORFS
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
scripts=/home/GLBRCORG/bpeterson26/HellsCanyon/code/generalUse
binsGood=~/HellsCanyon/dataEdited/binning/manualBinning/binsGood
binsGoodList=~/HellsCanyon/dataEdited/binning/manualBinning/goodBins_list.txt
cd $binsGood
mkdir ORFs

cat $binsGoodList | while read bin
do
  if [ ! -e ORFs/$bin.gff ]; then
    echo "Predicting proteins for" $bin
    prodigal -i DNA/$bin.fna \
              -o ORFs/$bin.gff \
              -f gff \
              -a ORFs/$bin.faa \
              -d ORFs/$bin.fna \
              -p single
  else
    echo $bin "already run."
  fi
done

# Clean up the gene sequence files
cd ORFs
cat $binsGoodList | while read bin
do
  echo "Cleaning" $bin
  python $scripts/cleanFASTA.py $bin.fna
  mv -f $bin.fna_temp.fasta $bin.fna
  python $scripts/cleanFASTA.py $bin.faa
  mv -f $bin.faa_temp.fasta $bin.faa
done

# Combine gene sequence files and generate S2B/G2B files
cd $binsGood
$scripts/Fasta_to_Scaffolds2Bin.sh -e fna \
                                    -i DNA \
                                    > binsGood_S2B.tsv
cat DNA/*.fna > DNA.fna
$scripts/Fasta_to_Scaffolds2Bin.sh -e faa \
                                    -i ORFs \
                                    > binsGood_G2B.tsv
cat ORFs/*.faa > ORFs.faa



####################################################
####################################################
# Run ANI comparisons on good bins
####################################################
####################################################
# The following workflow is the code to run
# Sarah Stevens's ANI calculator on a
# folder full of bins.
# Details: https://github.com/sstevens2/ani_compare_dag


mkdir ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/ANI_comparison
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/ANI_comparison
wget https://ani.jgi-psf.org/download_files/ANIcalculator_v1.tgz
tar -xzvf ANIcalculator_v1.tgz
git clone https://github.com/sstevens2/ani_compare_dag.git
mv ani_compare_dag HCC_bins_ANI
cd HCC_bins_ANI/
mkdir goodBins



cp ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/ORFs/*fna goodBins/
# from GLBRC to CHTC into ANI_comparison/ani_compare_dag/
echo 'goodBins' > groupslist.txt
# Change path of executable and
# transfer_input_files lines
#executable = /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binning/manualBinning/binsGood/ANI_comparison/HCC_bins_ANI/group.sh
#transfer_input_files = /home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_binning/binning_initial/binsGood/ANI_comparison/ANIcalculator_v1/ANIcalculator,/home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_binning/binning_initial/binsGood/ANI_comparison/ANIcalculator_v1/nsimscan,$(spllist),$(totransfer)
condor_submit_dag runAllANIcompare.dag
# Download output file (goodBins.all.ani.out.cleaned)
# to my computer:
# dataEdited/2019_binning/binning_initial/binsGood/goodBins.all.ani.out.cleaned

cd ~/Everglades/dataEdited/2019_binning/binning_initial
mkdir binsFinal
mkdir binsFinal/DNA
mkdir binsFinal/ORFs
cat binsFinal_list.txt | while read bin
do
  cp binsGood/DNA/$bin.fna binsFinal/DNA
  cp binsGood/ORFs/$bin* binsFinal/ORFs
done






####################################################
####################################################
# Metabolic analyses
####################################################
####################################################

##########################
# Custom set of metabolic HMMs
##########################

screen -S HCC_metabolic_HMMs
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate batch_HMMs
PYTHONPATH=''
PERL5LIB=''
binsGood=~/HellsCanyon/dataEdited/binning/manualBinning/binsGood
scripts=~/HellsCanyon/code/generalUse
metabolic_HMMs=~/HellsCanyon/references/metabolic_HMMs
cd $binsGood
mkdir metabolism
chmod +x $scripts/batch_HMMs.py

python $scripts/batch_HMMs.py --orf_file $binsGood/ORFs.faa \
                              --g2b $binsGood/binsGood_G2B.tsv \
                              --hmm_folder $metabolic_HMMs\
                              --hmm_csv $metabolic_HMMs.csv \
                              --output $binsGood/metabolism/batch_HMMs
conda deactivate


##########################
# Search for MHCs
##########################
mkdir $binsGood/metabolism/MHCs
cd $binsGood/metabolism/MHCs
$scripts/Find_multiheme_protein.py $binsGood/ORFs.faa 3
mv $binsGood/ORFs_3_heme* .

echo -e "binID\tgeneID\themeCount" > heme_count_bins.tsv
tail -n +2 ORFs_3_heme_count.txt | awk -F '\t' '{ print $1 }' | while read geneID
do
  binID=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' $binsGood/binsGood_G2B.tsv`
  hemeCount=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' ORFs_3_heme_count.txt`
  echo -e $binID"\t"$geneID"\t"$hemeCount
  echo -e $binID"\t"$geneID"\t"$hemeCount >> heme_count_bins.tsv
done

##########################
# Search for MHCs with 10 heme-binding sites
##########################
mkdir $binsGood/metabolism/MHCs_10
cd $binsGood/metabolism/MHCs_10
$scripts/Find_multiheme_protein.py $binsGood/ORFs.faa 10
mv $binsGood/ORFs_10_heme* .

echo -e "binID\tgeneID\themeCount" > heme_count_bins_10.tsv
tail -n +2 ORFs_10_heme_count.txt | awk -F '\t' '{ print $1 }' | while read geneID
do
  binID=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' $binsGood/binsGood_G2B.tsv`
  hemeCount=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' ORFs_10_heme_count.txt`
  echo -e $binID"\t"$geneID"\t"$hemeCount
  echo -e $binID"\t"$geneID"\t"$hemeCount >> heme_count_bins.tsv
done

##########################
# Search for BBOMPs
##########################
# Pull out names of adjacent genes
binsGood=~/HellsCanyon/dataEdited/binning/manualBinning/binsGood
cd $binsGood/metabolism
mkdir PCC
mkdir PCC/list
rm -f $binsGood/metabolism/PCC/adjacent_genes_all_list.txt

cat MHCs/ORFs_3_heme_list.txt | while read gene
do
  echo "Working on" $gene
  scaffold=$(echo $gene | rev | cut -d"_" -f2- | rev)
  assembly=$(echo $gene | rev | cut -d"_" -f3- | rev)
  ORFnumber=$(echo $gene | rev | cut -d"_" -f1 | rev)
  preceedingORFnumber=$(expr $ORFnumber - 1)
  followingORFnumber=$(expr $ORFnumber + 1)
  echo $scaffold"_"$preceedingORFnumber >> PCC/list/adjacent_genes_all_list.txt
  echo $scaffold"_"$followingORFnumber >> PCC/list/adjacent_genes_all_list.txt
done

# Find unique gene names
cd PCC
wc -l list/adjacent_genes_all_list.txt
sort list/adjacent_genes_all_list.txt | \
  uniq \
  > list/adjacent_genes_unique_list.txt

# Pull out adjacent genes
rm -f adjacent_genes.faa
cat list/adjacent_genes_unique_list.txt | while read geneID
do
  assembly=$(echo $geneID | rev | cut -d"_" -f3- | rev)
  echo "Looking for" $geneID "in" $assembly
  grep -A 1 -m 1 $geneID$ $binsGood/ORFs.faa >> adjacent_genes.faa
done


#########################
# Search adjacent genes for BBOMPs
#########################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""
cd $binsGood/metabolism/PCC

# Run the PCC HMM
pcc_omp_HMM=~/HellsCanyon/references/metabolicProteins/EET/pcc_omp.HMM
hmmsearch --tblout pcc_omp_custom.out \
          -T 40 \
          $pcc_omp_HMM \
          adjacent_genes.faa \
          > pcc_omp_custom_output.txt
# Pull out the sequences of interest
scripts=~/HellsCanyon/code/generalUse/
python $scripts/extract_protein_hitting_HMM.py \
        pcc_omp_custom.out \
        adjacent_genes.faa \
        pcc_omp_custom.faa


#########################
# Generate a phylogeny
#########################

cd $binsGood/metabolism/PCC
mkdir references

# Find references from RefSeq
blastp -query pcc_omp_custom.faa \
        -db ~/references/ncbi_db/refseq/refseq_protein \
        -evalue 0.005 \
        -outfmt '6 sseqid sseq' \
        -max_target_seqs 5 \
        -num_threads 3 \
        -out references/refseq_pcc_omp.txt
# Pull out amino acid sequences and dereplicate them
cd references
awk -F '\t' '{ print ">"$1"\n"$2 }' refseq_pcc_omp.txt > refseq_pcc_omp.faa
sed -i 's/ref|//' refseq_pcc_omp.faa
sed -i 's/|//' refseq_pcc_omp.faa
sed -i 's/-//g' refseq_pcc_omp.faa
# Dereplicate refseq sequences against reference dataset
cd-hit-2d -i ~/HellsCanyon/references/metabolicProteins/EET/pcc_omp.faa \
          -i2 refseq_pcc_omp.faa \
          -o blast_pcc_omp_uniq.faa \
          -c 0.97 \
          -n 5
# Get final list of refseq references
grep '>' blast_pcc_omp_uniq.faa | \
  sed 's/>//' \
  > blast_pcc_omp_uniq_list.txt
# Concatenate sequences
cd $binsGood/metabolism/PCC
mkdir phylogeny
cat references/blast_pcc_omp_uniq.faa \
    ~/HellsCanyon/references/metabolicProteins/EET/pcc_omp.faa \
    pcc_omp_custom.faa \
    > phylogeny/bbomp_phylogeny.faa

# Generate alignment
cd phylogeny
muscle -in bbomp_phylogeny.faa \
        -out bbomp_phylogeny.afa
# Trim the alignment
trimal -in bbomp_phylogeny.afa \
        -out bbomp_phylogeny_trimmed.afa \
        -gt 0.5
FastTree bbomp_phylogeny_trimmed.afa > bbomp.tree


# Retrieve reference information on local computer
HCC
cd dataEdited/binning/metabolism/PCC
epost -db protein -input blast_pcc_omp_uniq_list.txt | \
    esummary | \
    xtract -pattern DocumentSummary -element AccessionVersion,Organism > refseq_bbomp_metadata.tsv


#########################
# Check PCC_porin hits from batch HMMs
#########################
binsGood=~/HellsCanyon/dataEdited/binning/manualBinning/binsGood
cd $binsGood/metabolism
mkdir PCC_batchHMMs
mkdir PCC_batchHMMs/PCC

cp batch_HMMs/alignments/hmm_hits/pcc_porin.faa PCC_batchHMMs/PCC/pcc_porin.faa
cd PCC_batchHMMs/PCC/

grep '>' pcc_porin.faa | sed 's/>//' | awk '{ print $1 }' | while read gene
do
  #echo "Working on" $gene
  scaffold=$(echo $gene | rev | cut -d"_" -f2- | rev)
  assembly=$(echo $gene | rev | cut -d"_" -f3- | rev)
  ORFnumber=$(echo $gene | rev | cut -d"_" -f1 | rev)
  preceedingORFnumber=$(expr $ORFnumber - 1)
  followingORFnumber=$(expr $ORFnumber + 1)
  grep $scaffold"_"$preceedingORFnumber $binsGood/metabolism/MHCs/ORFs_3_heme_count.txt
  grep $scaffold"_"$followingORFnumber $binsGood/metabolism/MHCs/ORFs_3_heme_count.txt
done

#########################
# Confirm dsrA phylogeny
#########################
screen -S HCC_dsrA_tree
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=''
scripts=~/HellsCanyon/code/generalUse
binsGood=~/HellsCanyon/dataEdited/binning/manualBinning/binsGood

# Set up directory
cd $binsGood/metabolism
mkdir dsrA

# Copy over my sequences and align them
sed 's/*//' batch_HMMs/alignments/hmm_hits/dsrA.faa > dsrA/dsrA.faa
cd dsrA
grep '>' dsrA.faa | \
  sed 's/>//' \
  > dsrA_list.txt
muscle -in dsrA.faa \
        -out dsrA.afa

# Copy in reference sequences
cp ~/references/metabolicProteins/sulfur/dsrA/dsrA_karthik_clean.afa \
    dsrA_karthik_clean.afa

# Align seqs
muscle -profile \
        -in1 dsrA_karthik_clean.afa \
        -in2 dsrA.afa \
        -out dsrA_phylogeny.afa
# Trim alignment
trimal -in dsrA_phylogeny.afa \
        -out dsrA_phylogeny_trimmed.afa \
        -gt 0.5

# Generate ML tree
#screen -S EG_dsrA_tree
python $scripts/cleanFASTA.py dsrA_phylogeny_trimmed.afa
mv -f dsrA_phylogeny_trimmed.afa_temp.fasta dsrA_phylogeny_trimmed.afa
FastTree dsrA_phylogeny_trimmed.afa \
    > dsrA_phylogeny_trimmed.tree



#########################
# Confirm dsrA phylogeny
#########################
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/
mkdir metabolism/forMETABOLIC
ls DNA/*fna | sed 's/DNA\///' | sed 's/.fna//' | while read bin
do
  cp DNA/$bin.fna metabolism/forMETABOLIC/$bin.fasta
done



#########################
# Run METABOLIC
#########################
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/metabolism/forMETABOLIC
condor_submit ~/HellsCanyon/code/submission/METABOLIC_processing.sub



#########################
# Run Kofamscan on ORFs
#########################
cp -avr ~/references/kofamscan_DBs/kofam_scan-1.3.0 ~/HellsCanyon/code/
# Update config.yml file
binsGood=~/HellsCanyon/dataEdited/binning/manualBinning/binsGood
mkdir $binsGood/metabolism/KOFAM_output

# Get list of bins
cd $binsGood/ORFs
ls *fna | sed 's/.fna//' > $binsGood/bins_list.txt

# Run kofamscan on each set of bins
screen -S kofamscan
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate kofamscan
PYTHONPATH=""
PERL5LIB=""
binsGood=~/HellsCanyon/dataEdited/binning/manualBinning/binsGood
KOFAM_output=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binning/manualBinning/binsGood/metabolism/KOFAM_output
ORFs=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binning/manualBinning/binsGood/ORFs
cd /home/GLBRCORG/bpeterson26/HellsCanyon/code/kofam_scan-1.3.0

rm -f $KOFAM_output/output_forDecoder.tsv
cat $binsGood/bins_list.txt | while read binID
do
  if [ ! -d $KOFAM_output/$binID ]; then
    echo "Running KOFAMscan on" $binID
    ./exec_annotation -f detail-tsv \
                      -o $KOFAM_output/$binID.tsv \
                      $ORFs/$binID.faa

    echo "Cleaning up results from" $binID
    grep '*' $KOFAM_output/$binID.tsv > $KOFAM_output/$binID\_qualityHits.tsv
    binIDclean=`echo $binID | sed 's/_//g'`
    awk -F '\t' -v binIDclean="$binIDclean" '{ print binIDclean"_"$2"\t"$3 }' $KOFAM_output/$binID\_qualityHits.tsv >> $KOFAM_output/output_forDecoder.tsv
  fi
done
conda deactivate


#########################
# Run KEGG-decoder on output
#########################
screen -S kegg_decoder
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate kegg_decoder
PYTHONPATH=""
PERL5LIB=""
KOFAM_output=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binning/manualBinning/binsGood/metabolism/KOFAM_output
KEGG-decoder --input $KOFAM_output/output_forDecoder.tsv \
              --output $KOFAM_output/bin_metabolism_decoded_html.txt \
              --vizoption interactive


#########################
# Run genomes through Charles's workflow
#########################

cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/metabolism
mkdir FeGenie
cp /home/GLBRCORG/cnolmsted/FeGenie_EET_Stuff/PipelineOutput/All_EET_Proteins_bothMethods_HellsCanyon13genomes.csv \
    FeGenie


#########################
# Identify CAZymes in each bin
#########################
# On local computer
cd /Users/benjaminpeterson/Documents/research/HellsCanyon/dataEdited/binning/bin_orfs
split -l 60000 ORFs.faa ORFs_split_



#########################
# Characterize MoORs in bins
#########################
mkdir ~/HellsCanyon/references/custom_hmms
cp ~/5M/references/MoORs/MoOR.HMM ~/HellsCanyon/references/custom_hmms
cp ~/5M/references/MoORs/MoOR_final_reference.afa ~/HellsCanyon/references/custom_hmms/MoOR_reference.afa
# Search bins for MoORs
screen -S HCC_MoORs
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
scripts=~/HellsCanyon/code/generalUse/

mkdir metabolism/MoORs
hmmsearch --tblout metabolism/MoORs/MoOR.out \
          -T 50 \
          ~/HellsCanyon/references/custom_hmms/MoOR.HMM \
          ORFs.faa

# Pull out the gene names, use this to find correct portion of G2B file.
cd metabolism/MoORs/

rm -f MoOR_G2B.tsv
rm -f putative_MoORs.faa

grep -v '#' MoOR.out | \
  awk '{ print $1 }' | while read gene
  do
    awk -F '\t' -v gene="$gene" '$1 == gene { print $0 }' ../../binsGood_G2B.tsv >> MoOR_G2B.tsv
    grep -A 1 $gene\$ ../../ORFs.faa >> putative_MoORs.faa
    tail -n 1 MoOR_G2B.tsv
  done
sed -i 's/*//g' putative_MoORs.faa

hmmalign -o putative_MoORs.sto \
          ~/HellsCanyon/references/custom_hmms/MoOR.HMM \
          putative_MoORs.faa
$scripts/convert_stockhold_to_fasta.py putative_MoORs.sto
# Edit this in Geneious
muscle -profile \
        -in1 putative_MoORs_cut.afa \
        -in2 ~/HellsCanyon/references/custom_hmms/MoOR_reference.afa \
        -out putative_MoORs_ref_1.afa

python $scripts/cleanFASTA.py putative_MoORs_ref_1.afa
mv putative_MoORs_ref_1.afa_temp.fasta putative_MoORs_ref_1.afa
sed 's/-//g' putative_MoORs_ref_1.afa > putative_MoORs_ref_1.faa

# Mask the alignment at 50% gaps
trimal -in putative_MoORs_ref_1.afa \
        -out putative_MoORs_ref_1_masked.afa \
        -gt 0.5
FastTree putative_MoORs_ref_1_masked.afa > putative_MoORs_1.tree

raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s putative_MoORs_ref_1_masked.afa \
        -n hgcA


####################################################
####################################################
# hgcA analysis
####################################################
####################################################

screen -S HCC_hgcA_bins
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
binsGood=~/HellsCanyon/dataEdited/binning/manualBinning/binsGood
hgcaAssembly=~/HellsCanyon/dataEdited/binning/manualBinning/binsGood

scripts=~/HellsCanyon/code/generalUse/
cd $binsGood
mkdir hgcA

hmmsearch --tblout hgcA/hgcA.out \
          --cpu 4 \
          --cut_tc \
          ~/references/hgcA/hgcA.hmm \
          ORFs.faa \
          > hgcA/hgcA_report.txt
cd hgcA
grep -v "#" hgcA.out | awk '{ print $1 }' > hgcA_list.txt

# Get hgcA to bin file
echo -e 'hgcA_ID\tbinID' > hgcA_to_bin.tsv
cat hgcA_list.txt | while read hgcA
do
  awk -F '\t' -v hgcA="$hgcA" '$1 == hgcA { print $0 }' $binsGood/binsGood_G2B.tsv >> hgcA_to_bin.tsv
done



####################################################
####################################################
# Pull out coverage information from bam files
####################################################
####################################################

screen -S HCC_bin_coverage
binsGood=~/HellsCanyon/dataEdited/binning/manualBinning/binsGood
mkdir $binsGood/depth
cd /home/GLBRCORG/bpeterson26/HellsCanyon/code/
chmod +x executables/aggregate_depth_bins.sh
condor_submit submission/aggregate_depth_bins.sub
