#!/bin/sh

############################################
# code/prolixibacteraceae_details/gene_neighborhoods/prolixibacteraceae_GN.sh
# Benjamin D. Peterson

# This script will pull out the needed data
# to do the gene neighborhood analyses on
# hgcA.
############################################

if [ -f prolixibacteraceae_details ]; then
  rm -rf prolixibacteraceae_details
fi
working_directory=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/prolixibacteraceae_details
original_prolix_analysis=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny/prolixibacteraceae
HomeBio=/home/GLBRCORG/bpeterson26/HellsCanyon/code/homemade_bioinformatic_scripts
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh


############################################
# Get set up
############################################
echo "LOADING UP FOLDERS AND NEEDED GENOMES"
cd ~/HellsCanyon/dataEdited/binAnalysis

mkdir prolixibacteraceae_details
mkdir prolixibacteraceae_details/genomes

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
cd $working_directory
mkdir ORFs
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

cat genome_list.txt | while read accessionID
do
  bash $HomeBio/IMMA_ORF_STAN.sh -b $accessionID \
                                -i genomes/$accessionID.fna \
                                -o ORFs \
                                -c $HomeBio/fasta_stuff/cleanFASTA.py
done



############################################
# Pull out scaffolds of interest
############################################
echo "PULLING OUT GENE NEIGHBORHOOD INFORMATION FOR HGCA GENES"

conda deactivate
conda activate py_viz
PYTHONPATH=''
IFS=$'\n'
mkdir $working_directory/GN_of_hgcA
grep '>' $original_prolix_analysis/hgcA_hits/hgcA_raw.faa | sed 's/>//' > GN_of_hgcA/hgcA_id_list.txt

working_directory=/Users/benjaminpeterson/Documents/programs/homemade_bioinformatic_scripts/doctor_petersons_neighborhood/prolixibacteraceae_details
HomeBio=/Users/benjaminpeterson/Documents/programs/homemade_bioinformatic_scripts/

cd $working_directory
cat GN_of_hgcA/hgcA_id_list.txt | while read orfFastaID
do
  binID=`awk -F '\t' -v orfFastaID="$orfFastaID" '$1 == orfFastaID { print $2 }' ORFs_G2B.tsv`
  python $HomeBio/doctor_petersons_neighborhood.py --bin_id $binID \
                                                   --bin_file genomes/$binID.fna \
                                                   --orf_fasta_id $orfFastaID \
                                                   --ORF_location ORFs \
                                                   --size_of_block 5000 \
                                                   --outputLocation GN_of_hgcA
done




cd GN_of_hgcA
# rm -f hgcA_geneNeighborhood.gff hgcA_geneNeighborhood.fna
cat *_neighborhood.gff > hgcA_geneNeighborhood_all.gff
cat *_neighborhood.fna > hgcA_geneNeighborhood_all.fna
rm -f *_neighborhood.*
