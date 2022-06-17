#!/bin/sh

##################################################
##################################################
# code/binning/PVC_details/metabolic_genes_PVC.sh
# Benjamin D. Peterson

# This script will look for MHC content and metabolic
# genes that we include in the batch HMM analysis
# to see if there are major metabolic differences
# between the PVC families, particularly the two
# Kiritimatiellaeota families.
##################################################
##################################################


working_directory=~/HellsCanyon/dataEdited/binAnalysis/PVC_details
HomeBio=/home/GLBRCORG/bpeterson26/HellsCanyon/code/HomeBio
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate gtdbtk
PYTHONPATH=''
PERL5LIB=''

cd $working_directory
mkdir taxonomy



####################################################
####################################################
# Run genomes through GTDB
####################################################
####################################################
gtdbtk classify_wf \
        --cpus 16 \
        --extension fna \
        --genome_dir ./ORFs \
        --prefix PVC_taxonomy_GTDB \
        --out_dir taxonomy
# Summarize them
cd taxonomy
awk -F '\t' '{ print $1"\t"$2 }' gtdbtk.*.summary.tsv \
        > taxonomy_summary.txt
conda deactivate
