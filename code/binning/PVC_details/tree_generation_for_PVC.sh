#!/bin/sh

##################################################
##################################################
# code/binning/phylogenies/PVC_tree.sh
# Benjamin D. Peterson
##################################################
##################################################


##################################################
# Set up paths
##################################################
PVC_details=~/HellsCanyon/dataEdited/binAnalysis/PVC_details
HomeBio=/home/GLBRCORG/bpeterson26/HellsCanyon/code/HomeBio
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

mkdir $PVC_details/tree_generation
cd $PVC_details

##################################################
# Run tree generation workflow
##################################################
python $HomeBio/bin/gene_alignment_from_bins.py \
          --analysis_name HCC_PVC \
          --orf_file HCC_PVC_ORFs.faa \
          --g2b_file HCC_PVC_ORFs_G2B.tsv \
          --hmm_list $HomeBio/databases/HMMs/rp16_bact_HMM_list.txt \
          --hmm_location $HomeBio/databases/HMMs/hmm_folder \
          --output_location tree_generation \
          --threads_to_use 16 \
          --minimum_hits 8 \
          --masking_threshold 0.5 \
          --tree_program RAxML
