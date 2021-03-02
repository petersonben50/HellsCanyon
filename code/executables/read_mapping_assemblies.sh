#!/bin/sh


##############################
# code/executables/read_mapping_assemblies.sh
# Benjamin D. Peterson

# This script uses bowtie2 to map reads to the
# assemblies for 2018 analyses.
# It also cleans and processes the output using
# samtools

# Submission file loops over all the assemblies.
##############################


# Will read in these variables from
# submission file:
# $assembly
# $metagenome

# Will read in these paths from the submission file:
# read_storage: This is the folder where the cleaned reads are located.
# scaffolds: This is the folder where the scaffolds that will be mapped to are located
# output: This is the folder where the output mapping folder will go. This script cds into that folder.
# metagenomeList: This is a list of the metagenomes that we want mapped to each assembly


##############################
# Set up the environment
##############################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PERL5LIB=""
cd $output



##############################
# Run the mapping
##############################

if [ -e indices/bowtie_index_$assembly.1.bt2 ]; then
  if [ ! -e $metagenome\_to_$assembly.bam ]; then
    echo "Mapping" $metagenome "to" $assembly
    bowtie2 -x indices/bowtie_index_$assembly \
            -1 $read_storage/$metagenome\_R1.fastq.gz \
            -2 $read_storage/$metagenome\_R2.fastq.gz \
            -U $read_storage/$metagenome\_merged.fastq.gz,$read_storage/$metagenome\_single.fastq.gz \
            -p 12 \
            -S $metagenome\_to_$assembly.sam
    samtools view $metagenome\_to_$assembly.sam \
                  -o $metagenome\_to_$assembly.unsorted.bam
    # Sort the BAM files
    samtools sort -m 10G \
                  -@ 8 \
                  $metagenome\_to_$assembly.unsorted.bam \
                  -o $metagenome\_to_$assembly.bam
    # Index the BAM files.
    samtools index $metagenome\_to_$assembly.bam

    # Clean up
    rm -f $metagenome\_to_$assembly.unsorted.bam
    rm -f $metagenome\_to_$assembly.sam
  else
    echo "Mapping of" $metagenome "to" $assembly "is already done"
  fi
else
  echo "Need to complete indexing on this assembly:" $assembly "to map" $metagenome
fi
