#!/bin/sh

##############################
# code/executables/nonpareil3_run.sh
# Benjamin D. Peterson

# This is a general script to run nonpareil,
# given a single metagenome (no paired reads).
##############################


#########################
# Set up
#########################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
conda activate nonpareil3
cd $workingDirectory

#########################
# Copy over and unzip fastq file
#########################
cp $metagenomeFolder/$metagenome\_R1.fastq.gz ./
gunzip $metagenome\_R1.fastq.gz


#########################
# Run Nonpareil
#########################
nonpareil -s $metagenome\_R1.fastq \
          -T kmer \
          -f fastq \
          -t 8 \
          -X 100000 \
          -k 32 \
          -b $metagenome\_NPoutput


#########################
# Clean up
#########################
rm -f $metagenome\_R1.fastq
