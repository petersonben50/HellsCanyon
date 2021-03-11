#!/bin/sh

##################################################
##################################################
# code/executables/aggregate_depth_proteins.sh
# Benjamin D. Peterson

# This script will calculate the depth of coverage
# of a given set of scaffolds within a given set
# of metagenomes.
##################################################
##################################################

# Variables needed:
# metagenome
# depthDirectory
# mappingFolder
# scripts

#########################
# Set up for samtools calculation of depth
#########################

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
rm -f $depthDirectory/$metagenome\_depth_raw.tsv
conda activate bioinformatics
PERL5LIB=""


#########################
# Calculate depth over each residue in each scaffold
#########################
cat $depthDirectory/scaffold_all_list.txt | while read scaffold
do
  assembly=$(echo $scaffold | awk -F '_' '{ print $1 }')
  if [ -e $mappingFolder/$metagenome\_to_$assembly.bam ]; then
    echo "Calculating coverage of" $metagenome "over" $scaffold
    samtools depth -a -r $scaffold $mappingFolder/$metagenome\_to_$assembly.bam \
        >> $depthDirectory/$metagenome\_depth_raw.tsv
  else
    echo $metagenome "was not mapped to" $scaffold
  fi
done

conda deactivate

#########################
# Average the depth of each residue over the entire
# scaffold, minus the first and last 150 base pairs
#########################
echo "Aggregating" $scaffold "depth information for" $metagenome
conda activate py_viz
PYTHONPATH=""
python $scripts/calculate_depth_contigs.py \
          $depthDirectory/$metagenome\_depth_raw.tsv \
          150 \
          $depthDirectory/$metagenome\_depth.tsv

rm -f $depthDirectory/$metagenome\_depth_raw.tsv
