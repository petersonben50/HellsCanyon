#!/bin/sh

######################
# code/metabolic_analyses/dsrA_analysis.sh
# Benjamin D. Peterson

# This set of scripts will contain the bash
# code we need to pull out the dsrA sequences
# from our assembly and process them.
######################


############################################
############################################
# Identify dsrA sequences with GID
############################################
############################################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PATH="/home/GLBRCORG/bpeterson26/BLiMMP/code/HomeBio/bin:$PATH"

#cd ~/BLiMMP/code
#git clone https://github.com/petersonben50/HomeBio
#cd ~/BLiMMP/code/HomeBio
#git pull
HomeBio=~/BLiMMP/code/HomeBio

# Set up dataset
# rm -fr ~/HellsCanyon/dataEdited/metabolic_analyses/dsrA_GID
python3 $HomeBio/bin/ABA_GID.py --orf_folder ~/HellsCanyon/dataEdited/assemblies/ORFs \
                  --hmm $HomeBio/reference_data/HMMs/hmm_folder/TIGR02064.HMM \
                  --output_location /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/metabolic_analyses/dsrA_GID \
                  --output_prefix dsrA \
                  --metagenome_list /home/GLBRCORG/bpeterson26/HellsCanyon/metadata/lists/MG_list_all.txt \
                  --metagenomes_location /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/mapping \
                  --reference_aa_dataset $HomeBio/reference_data/sequence_databases/dsrA/muller_DsrA_dataset_final.faa \
                  --number_threads 30 \
                  --cluster_cutoff 0.8 \
                  > ~/HellsCanyon/dataEdited/metabolic_analyses/GID_log_dsrA.txt

