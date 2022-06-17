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
metabolic_HMMs=~/HellsCanyon/references/metabolic_HMMs
scripts=~/HellsCanyon/code/generalUse
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

cd $working_directory
mkdir metabolism


####################################################
####################################################
# Custom set of metabolic HMMs
####################################################
####################################################

#chmod +x $scripts/batch_HMMs.py
python $HomeBio/batch_HMMs/batch_HMMs.py --orf_file $working_directory/HCC_PVC_ORFs.faa \
                                          --g2b $working_directory/HCC_PVC_ORFs_G2B.tsv \
                                          --hmm_folder $metabolic_HMMs\
                                          --hmm_csv $metabolic_HMMs.csv \
                                          --output $working_directory/metabolism/batch_HMMs
conda deactivate



####################################################
####################################################
# Search for MHCs
####################################################
####################################################

mkdir $working_directory/metabolism/MHCs
cd $working_directory/metabolism/MHCs
$scripts/Find_multiheme_protein.py $working_directory/HCC_PVC_ORFs.faa 3
mv $working_directory/HCC_PVC_ORFs_3_heme* .

echo -e "binID\tgeneID\themeCount" > heme_count_bins.tsv
tail -n +2 HCC_PVC_ORFs_3_heme_count.txt | awk -F '\t' '{ print $1 }' | while read geneID
do
  binID=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' $working_directory/HCC_PVC_ORFs_G2B.tsv`
  hemeCount=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' HCC_PVC_ORFs_3_heme_count.txt`
  echo -e $binID"\t"$geneID"\t"$hemeCount
  echo -e $binID"\t"$geneID"\t"$hemeCount >> heme_count_bins.tsv
done
