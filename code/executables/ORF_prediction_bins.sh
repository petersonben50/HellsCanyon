#!/bin/sh


##############################
# code/executables/ORF_prediction_bins.sh
# Benjamin D. Peterson

# This script uses Prodigal to predict
# the open reading frames of the assemblies
# from the Hells Canyon metagenomes.
# Submission file loops over all the assemblies.
##############################


##############################
# Set up the environment
##############################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
cd /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binning/autoBinning/hqBinSet





##############################
# Run ORF prediction
##############################
if [ ! -e ORFs/$binID.gff ]; then
  echo "Predicting proteins for" $binID
  prodigal -i DNA/$binID.fna \
            -o ORFs/$binID.gff \
            -f gff \
            -a ORFs/$binID.faa \
            -d ORFs/$binID.fna \
            -p single
  python /home/GLBRCORG/bpeterson26/HellsCanyon/code/generalUse/cleanFASTA.py ORFs/$binID.fna
  mv -f ORFs/$binID.fna_temp.fasta ORFs/$binID.fna
  python /home/GLBRCORG/bpeterson26/HellsCanyon/code/generalUse/cleanFASTA.py ORFs/$binID.faa
  mv -f ORFs/$binID.faa_temp.fasta ORFs/$binID.faa
else
  echo "Already predicted ORFs for" $binID
fi
