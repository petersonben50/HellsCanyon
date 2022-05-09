#!/bin/sh

##################################################
##################################################
# code/binning/PVC_details/hgcA_analysis_for_PVC.sh
# Benjamin D. Peterson
##################################################
##################################################

working_directory=~/HellsCanyon/dataEdited/binAnalysis/PVC_details
HomeBio=/home/GLBRCORG/bpeterson26/HellsCanyon/code/HomeBio
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

mkdir $working_directory/hgcA_analysis
mkdir $working_directory/hgcA_analysis/hgcA_hits

##################################################
# Identify hgcA sequences
##################################################
cd $working_directory
python $HomeBio/bin/protein_identification_and_alignment.py \
    --analysis_name PVC_hgcA \
    --hmm_file ~/references/hgcA/hgcA.hmm \
    --orf_file $working_directory/HCC_PVC_ORFs.faa \
    --output_location $working_directory/hgcA_analysis/hgcA_hits \
    --cpus_to_use 10 \
    --hmm_cutoff TC

# After checking alignment to ensure they all belong.
# Make list
grep '>' $working_directory/hgcA_analysis/hgcA_hits/PVC_hgcA.afa | \
  sed 's/>//' > $working_directory/hgcA_analysis/PVC_hgcA_geneID_list.txt
# Generate G2B file for gene of interest
if [ -e $working_directory/hgcA_analysis/PVC_hgcA_G2B.tsv ]; then
  rm -f $working_directory/hgcA_analysis/PVC_hgcA_G2B.tsv
fi
cat $working_directory/hgcA_analysis/PVC_hgcA_geneID_list.txt | while read geneID; do
  awk -v geneID="$geneID" -F '\t' '$1 == geneID { print $0 }' $working_directory/HCC_PVC_ORFs_G2B.tsv \
      >> $working_directory/hgcA_analysis/PVC_hgcA_G2B.tsv
done


##################################################
# Generate hgcA tree
##################################################
mkdir $working_directory/hgcA_analysis/tree_generation
trimal -in $working_directory/hgcA_analysis/hgcA_hits/PVC_hgcA.afa \
       -out $working_directory/hgcA_analysis/tree_generation/PVC_hgcA_masked.afa \
       -gt 0.5
raxmlHPC-PTHREADS -f a \
                  -p 54457 \
                  -m PROTGAMMAAUTO \
                  -N autoMRE \
                  -x 2381 \
                  -T 16 \
                  -w $working_directory/hgcA_analysis/tree_generation \
                  -s $working_directory/hgcA_analysis/tree_generation/PVC_hgcA_masked.afa \
                  -n PVC_hgcA



##################################################
# Set up gene neighborhood analysis
##################################################
