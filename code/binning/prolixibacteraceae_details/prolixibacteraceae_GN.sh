#!/bin/sh

############################################
# code/prolixibacteraceae_details/gene_neighborhoods/prolixibacteraceae_GN.sh
# Benjamin D. Peterson

# This script will pull out the needed data
# to do the gene neighborhood analyses on
# hgcA.
############################################

############################################
# Prepare variables, directories, conda
############################################
cd /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis
if [ -f prolixibacteraceae_details ]; then
  rm -rf prolixibacteraceae_details
fi
mkdir prolixibacteraceae_details

working_directory=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/prolixibacteraceae_details
original_prolix_analysis=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny/prolixibacteraceae
HomeBio=/home/GLBRCORG/bpeterson26/HellsCanyon/code/HomeBio
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh

cd $working_directory
mkdir reports/



############################################
# Retrieve needed genomes
############################################
echo "LOADING UP FOLDERS AND NEEDED GENOMES"
mkdir genomes

grep '>' $original_prolix_analysis/hgcA_hits/hgcA_raw.faa | sed 's/>//' | while read geneID
do
  binID=`awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $2 }' $original_prolix_analysis/ORFs_G2B.tsv`
  cp $original_prolix_analysis/reference_bins/$binID.fna prolixibacteraceae_details/genomes
done
cp ~/HellsCanyon/dataEdited/binning/bins_hgcA_keepers/DNA/anvio_hgcA_0130.fna prolixibacteraceae_details/genomes
cd prolixibacteraceae_details/genomes
ls *fna | sed 's/.fna//' > ../genome_list.txt




############################################
# Predict ORFs
############################################
echo "PREDICTING ORFS WITH IMMA_ORF_STAN"
mkdir ORFs
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

# Predict ORFs
cat genome_list.txt | while read accessionID
do
  bash $HomeBio/PM4_binGeneration/IMMA_ORF_STAN_BINS.sh -b $accessionID \
                                                        -i genomes/$accessionID.fna \
                                                        -o ORFs \
                                                        -c $HomeBio/fasta_manipulation/cleanFASTA.py
done > reports/IMMA_ORF_STAN_log.txt

# Generate scaffold to bin file
bash $HomeBio/fasta_manipulation/Fasta_to_Scaffolds2Bin.sh -e faa \
                                                           -i $working_directory/ORFs \
                                                           > $working_directory/ORFs_G2B.tsv
cat $working_directory/ORFs/*.faa > $working_directory/ORFs.faa
conda deactivate




############################################
# Pull out scaffolds of interest
############################################
echo "PULLING OUT GENE NEIGHBORHOOD INFORMATION FOR HGCA GENES"

conda activate py_viz
PYTHONPATH=''
IFS=$'\n'
mkdir $working_directory/GN_of_hgcA
grep '>' $original_prolix_analysis/hgcA_hits/hgcA_raw.faa | sed 's/>//' > $working_directory/GN_of_hgcA/hgcA_id_list.txt

cat GN_of_hgcA/hgcA_id_list.txt | while read orfFastaID
do
  binID=`awk -F '\t' -v orfFastaID="$orfFastaID" '$1 == orfFastaID { print $2 }' ORFs_G2B.tsv`
  python $HomeBio/AM3_binBasedAnalyses/doctor_petersons_neighborhood.py --bin_id $binID \
                                                                         --bin_file genomes/$binID.fna \
                                                                         --orf_fasta_id $orfFastaID \
                                                                         --ORF_location ORFs \
                                                                         --size_of_block 5000 \
                                                                         --outputLocation GN_of_hgcA
done > reports/doctor_petersons_neighborhood_log.txt

# Concatenate the output files
cd GN_of_hgcA
rm -f hgcA_geneNeighborhood_all_orfs.faa
cat *_neighborhood.gff > hgcA_geneNeighborhood_all.gff
cat *_neighborhood.fna > hgcA_geneNeighborhood_all.fna
cat *.faa > hgcA_geneNeighborhood_all_orfs.faa
rm -f *_[0-9]*.faa
rm -f *_neighborhood.*
conda deactivate


############################################
# Run KOFAMSCAN on the neighborhood ORFs
############################################
echo "PULLING OUT GENE NEIGHBORHOOD INFORMATION FOR HGCA GENES"
conda activate kofamscan
PYTHONPATH=""
PERL5LIB=""
analysis_name=prolix_hgcA_GN
input_orfs=GN_of_hgcA/hgcA_geneNeighborhood_all_orfs.faa
output_location=GN_of_hgcA/kofamscan_output
cd $working_directory

rm -f $output_location/$analysis_name\_raw.tsv
echo "Running KOFAMscan for" $analysis_name
~/references/kofamscan_files/kofam_scan-1.3.0/exec_annotation -f detail-tsv \
                                                              -o $output_location/$analysis_name\_raw.tsv \
                                                              $input_orfs



############################################
# Cluster ORFs with CD-HIT
############################################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""

cd $working_directory
mkdir GN_of_hgcA/clustering
cd-hit -g 1 \
        -i GN_of_hgcA/hgcA_geneNeighborhood_all_orfs.faa \
        -o GN_of_hgcA/clustering/hgcA_geneNeighborhood_all_orfs_cluster_40.faa \
        -c 0.40 \
        -n 2 \
        -d 0
clstr2txt.pl GN_of_hgcA/clustering/hgcA_geneNeighborhood_all_orfs_cluster_40.faa.clstr \
  > GN_of_hgcA/clustering/hgcA_geneNeighborhood_all_orfs_cluster_40.tsv


# Files to download:
# GN_of_hgcA/clustering/hgcA_geneNeighborhood_all_orfs_cluster_40.tsv
# $output_location/$analysis_name\_quality_cleaned.tsv
# $working_directory/GN_of_hgcA/hgcA_id_list.txt



############################################
# Generate figure
############################################
"""
# On local machine
conda activate py_viz
PYTHONPATH=""
PERL5LIB=""

working_directory=/Users/benjaminpeterson/Documents/research/HellsCanyon/dataEdited/bins/binAnalysis/prolixibacteraceae_details
HomeBio=/Users/benjaminpeterson/Documents/programs/HomeBio

GFF_FILE = '/Users/benjaminpeterson/Documents/research/HellsCanyon/dataEdited/bins/binAnalysis/prolixibacteraceae_details/GN/hgcA_geneNeighborhood_all.gff'
ORF_DATA = '/Users/benjaminpeterson/Documents/research/HellsCanyon/dataEdited/bins/binAnalysis/prolixibacteraceae_details/GN/hgcA_geneNeighborhood_info.txt'

python $HomeBio/AM3_binBasedAnalyses/doctor_petersons_neighborhood_visualization.py --gff_file $working_directory/GN/hgcA_geneNeighborhood_all.gff \
                                                                                    --orf_data $working_directory/GN/hgcA_geneNeighborhood_info.txt \
                                                                                    --output_location $working_directory/GN/hgcA_geneNeighborhood_plot.pdf

"""
