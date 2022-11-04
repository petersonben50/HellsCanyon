#!/bin/sh

##############################
# code/executables/emirge_run.sh
# Benjamin D. Peterson

# This is a general script to run EMIRGE
# on a specified metagenome.
##############################


#########################
# Set up
#########################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate emirge
PYTHONPATH=''
cd $emirge_output
PATH=$PATH:/opt/bifxapps/bin/

#########################
# Unzip fastq file
#########################
gunzip < $read_storage/$metagenome\_R1.fastq.gz \
    > $emirge_output/$metagenome\_R1.fastq
gunzip < $read_storage/$metagenome\_R2.fastq.gz \
    > $emirge_output/$metagenome\_R2.fastq

#########################
# Set parameters
#########################
insert_size=$(awk -F '\t' -v metagenome="$metagenome" '$1 == metagenome { print $2}' /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/16S_from_MG/metagenome_inserts.txt)
insert_sd=$(awk -F '\t' -v metagenome="$metagenome" '$1 == metagenome { print $4}' /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/16S_from_MG/metagenome_inserts.txt)

#########################
# Run emirge
#########################
echo "Running EMIRGE on" $metagenome
emirge.py $emirge_output/$metagenome \
          -1 $emirge_output/$metagenome\_R1.fastq \
          -2 $emirge_output/$metagenome\_R2.fastq \
          -f /home/GLBRCORG/bpeterson26/references/emirge_db/SILVA_138_SSURef_NR99_tax_silva_trunc.ge1200bp.le2000bp.0.97.fixed.fasta \
          -b /home/GLBRCORG/bpeterson26/references/emirge_db/SILVA_138_SSURef_NR99_tax_silva_trunc.ge1200bp.le2000bp.0.97.fixed \
          -l 150 \
          -i $insert_size \
          -s $insert_sd \
          -a 8 \
          --phred33

#########################
# Clean up
#########################
rm -f $emirge_output/$metagenome\_R1.fastq \
      $emirge_output/$metagenome\_R2.fastq
