#!/bin/sh

#########################
# executables/automatic_binning.sh
# Benjamin D. Peterson

# This script
#########################


##########################
# Activate conda environment
##########################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate binning
PYTHONPATH=''
PERL5LIB=''
R_LIBS=''


####################################################
####################################################
# Generate bins from Metabat2
####################################################
####################################################

##########################
# Get set up
##########################
cd $metabatOutput
mkdir $assembly\_output

##########################
# Generate list of metagenomes to include in binning
##########################
awk -F '\t' -v assembly="$assembly" '$2 == assembly { print $1 }' $mappingKey > metagenomes_for_$assembly\_mapping.txt

##########################
# Generate bins from metabat2
##########################
# First need to generate a depth file
if [ ! -e $assembly\_output/depth_to_$assembly.txt ]; then
  echo "Summarizing depths for" $assembly
  python $scripts/depth_file_generation_metabat.py metagenomes_for_$assembly\_mapping.txt \
                                                    $mappingLocation \
                                                    $assembly \
                                                    $assembly\_output
else
  echo "Depth file generation done for" $assembly
fi
# Then run the binning
if [ ! -e $assembly\_metabat* ]; then
  echo "Binning" $assembly
  metabat2 -i $scaffoldsLocation/$assembly\_filtered_scaffolds.fna \
            -a $assembly\_output/depth_to_$assembly.txt \
            -o $assembly\_output/$assembly\_metabat \
            -m 2000
else
  echo "Binning of" $assembly "already done"
fi

##########################
# Rename bins from metabat2
##########################
cd $assembly\_output
ls *fa | while read bin
do
  newBinName=$(echo $bin | \
                sed "s/$assembly\_metabat./$assembly\_metabat_/" | \
                sed "s/.fa/.fna/")
  echo "Renaming" $bin "to" $newBinName
  mv $bin $newBinName
done
echo "Generating S2B file for" $assembly
$scripts/Fasta_to_Scaffolds2Bin.sh -e fna > ../$assembly\_metabat_S2B.tsv
rm -f metagenomes_for_$assembly\_mapping.txt


####################################################
####################################################
# Generate bins using MaxBin2
####################################################
####################################################

##########################
# Get set up
##########################
cd $maxbinOutput
mkdir $assembly\_output

##########################
# Generate list of metagenomes to include in binning
##########################
awk -F '\t' -v assembly="$assembly" '$2 == assembly { print $1 }' $mappingKey > metagenomes_for_$assembly\_mapping.txt

##########################
# Prep depth files
##########################
# Generate abundance file from BAM files
cat metagenomes_for_$assembly\_mapping.txt | while read metagenome
do
  echo "Calculating depth of" $metagenome "over" $assembly
  samtools depth -aa $mappingLocation/$metagenome\_to_$assembly.bam \
    > $assembly\_output/$metagenome\_to_$assembly\_raw.depth
  python $scripts/calculate_depth_contigs.py $assembly\_output/$metagenome\_to_$assembly\_raw.depth \
                                              150 \
                                              $assembly\_output/$metagenome\_to_$assembly.depth
  rm -f $assembly\_output/$metagenome\_to_$assembly\_raw.depth
done
cd $assembly\_output/
ls *_to_$assembly.depth > abund_list.txt

#################################
# Create list of abundance files
#################################
run_MaxBin.pl -contig $scaffoldsLocation/$assembly\_filtered_scaffolds.fna \
              -abund_list abund_list.txt \
              -out $assembly\_maxbin \
              -thread 10
rm -f metagenomes_for_$assembly\_mapping.txt


##########################
# Rename bins from MaxBin2
##########################
ls *fasta | while read bin
do
  newBinName=$(echo $bin | \
                sed "s/$assembly\_maxbin./$assembly\_maxbin_/" | \
                sed "s/fasta/fna/")
  echo "Renaming" $bin "to" $newBinName
  mv $bin $newBinName
done
echo "Generating S2B file for" $assembly
$scripts/Fasta_to_Scaffolds2Bin.sh -e fna > ../$assembly\_maxbin_S2B.tsv


####################################################
####################################################
# Aggregate bins using Das Tool
####################################################
####################################################
cd $dasToolOutput
mkdir $assembly\_output
DAS_Tool -i $metabatOutput/$assembly\_metabat_S2B.tsv,$maxbinOutput/$assembly\_maxbin_S2B.tsv \
        -l metabat,maxbin \
        -c $scaffoldsLocation/$assembly\_filtered_scaffolds.fna \
        -o $assembly\_output/$assembly\_bins \
        --threads 9 \
        --write_bins 1 \
        --search_engine diamond \
        --create_plots 1 \
        --score_threshold 0.4
