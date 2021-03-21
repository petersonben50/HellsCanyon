#!/bin/sh

#########################
# executables/anvio_profiling.sh
# Benjamin D. Peterson

# This script takes the reads from each
# metagenome that have been mapped to the
# assemblies from this analysis and
# profiles them into an anvi'o database.
#########################

# Will read in these variables from
# submission file:
# $assembly

# Will read in these paths from submission file:
# output: Output location. Should be same place as the anvio databases
# mappingKey: tsv file that links metagenomes to assemblies
# mappingLocation: Location of the bam files from the mapping


##########################
# Activate conda environment
##########################
cd $output
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""


#########################
# Set up list of metagenomes for mapping
#########################
awk -F '\t' -v assembly="$assembly" '$2 == assembly { print $1 }' $mappingKey > metagenomes_for_$assembly\_mapping.txt


#########################
# Generate read profiles
#########################
cat metagenomes_for_$assembly\_mapping.txt | while read metagenome
do
  if [ ! -d $metagenome\_to_$assembly.profile ]; then
    echo "Profiling reads mapped from:" $metagenome "to" $assembly
    anvi-profile -c $assembly.db \
                  -i $mappingLocation/$metagenome\_to_$assembly.bam \
                  --write-buffer-size 500 \
                  --min-contig-length 2000 \
                  --num-threads 10 \
                  -o $metagenome\_to_$assembly.profile
  else
    echo "STOP: Profiling reads mapped to:" $assembly "from" $metagenome "is already done"
  fi
done


#########################
# Merge read profiles
#########################
if [ ! -d $assembly.merged ]; then
  anvi-merge *to_$assembly.profile/PROFILE.db \
              -o $assembly.merged \
              -c $assembly.db \
              -S $assembly\_merged \
              --skip-hierarchical-clustering
else
  echo "Merging profiles complete for" $assembly
fi


#########################
# Clean up
#########################
rm -f metagenomes_for_$assembly\_mapping.txt
