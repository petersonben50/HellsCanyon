#!/bin/sh

#########################
# executables/anvio_DB_prep.sh
# Benjamin D. Peterson

# This script will use anvi'o v6.2
# to generate anvi'o databases for
# a given assembly and will populate
# it with HMMs for SCG genes and
# hgcA.
#########################

##########################
# Activate conda environment
##########################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""

if [ ! -f $outputLocation/$assembly.db ]; then

  echo "Generating contig database for" $assembly

  ##########################
  # Generate anvio databases
  ##########################
  anvi-gen-contigs-database -f $scaffoldsLocation/$assembly\_filtered_scaffolds.fna \
                            -o $outputLocation/$assembly.db \
                            -n $assembly


  ##########################
  # Populate anvi'o database with HMMs
  ##########################

  # Run default anvio HMMs
  anvi-run-hmms -c $outputLocation/$assembly.db \
                --num-threads 7

  # Run set of custom HMMs
  anvi-run-hmms -c $outputLocation/$assembly.db \
                -H $customHMMs \
                --num-threads 7

  # Add taxonomic information
  anvi-run-scg-taxonomy -c $outputLocation/$assembly.db

else
  echo "Contig database for" $assembly "already exists."
fi
