#!/bin/sh

############################################
# code/prolixibacteraceae_details/prolixibacteraceae_GN.sh
# Benjamin D. Peterson

# This script will pull out the needed data
# to do the gene neighborhood analyses on
# hgcA.
############################################

############################################
# Prepare variables, directories, conda
############################################
cd /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis
#mv prolixibacteraceae_details prolixibacteraceae_details_take1
if [ -f prolixibacteraceae_details ]; then
  rm -rf prolixibacteraceae_details
fi
mkdir prolixibacteraceae_details

working_directory=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/prolixibacteraceae_details
original_prolix_analysis=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny/prolixibacteraceae
original_bact_analysis=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny/Bacteroidetes
HomeBio=/home/GLBRCORG/bpeterson26/HellsCanyon/code/HomeBio
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh

cd $working_directory
mkdir reports/



############################################
# Retrieve needed genomes
############################################
echo "LOADING UP FOLDERS AND NEEDED GENOMES"
mkdir genomes

grep '>' $original_prolix_analysis/hgcA_hits/hgcA_raw.afa | sed 's/>//' | while read geneID
do
  binID=`awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $2 }' $original_prolix_analysis/ORFs_G2B.tsv`
  if [ $binID != "anvio_hgcA_0130" ]; then
    echo "Moving" $binID
    cp $original_prolix_analysis/reference_bins/$binID.fna genomes
  fi
done

# Copy hgcA+ Bacteroidetes genomes
grep '>' $original_bact_analysis/hgcA_hits/hgcA_raw.faa  | sed 's/>//' | while read geneID
do
  binID=`awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $2 }' $original_bact_analysis/ORFs_G2B.tsv`
  if [ $binID == "anvio_hgcA_0130" ]; then
    echo "Moving" $binID "from this study"
    cp ~/HellsCanyon/dataEdited/binning/bins_hgcA_keepers/DNA/$binID.fna genomes
  elif echo $binID | grep -q "BAC_" ; then
    echo "Moving" $binID "from Mendota study"
    cp ~/5M/dataEdited/binAnalysis/bins_processed/$binID.fna genomes
  else
    echo "Moving" $binID "from references"
    cp $original_bact_analysis/scaffolds/$binID.fna genomes
  fi
done
cd genomes
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
# Look for genes of interest
############################################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
mkdir hgcA_hits

hmmsearch --tblout $working_directory/hgcA_hits/hgcA.out \
          --cpu 4 \
          --cut_tc \
          ~/references/hgcA/hgcA.hmm \
          $working_directory/ORFs.faa \
          > $working_directory/hgcA_hits/hgcA_HMM_output.txt
python $HomeBio/fasta_manipulation/extract_protein_hitting_HMM.py \
                $working_directory/hgcA_hits/hgcA.out \
                $working_directory/ORFs.faa \
                $working_directory/hgcA_hits/hgcA_raw.faa
grep '>' $working_directory/hgcA_hits/hgcA_raw.faa | \
  sed 's/>//' \
  > hgcA_hits/hgcA_id_list.txt
conda deactivate


############################################
# Pull out scaffolds of interest
############################################
echo "PULLING OUT GENE NEIGHBORHOOD INFORMATION FOR HGCA GENES"

conda activate py_viz
PYTHONPATH=''
IFS=$'\n'
mkdir $working_directory/GN_of_hgcA


cat hgcA_hits/hgcA_id_list.txt | while read orfFastaID
do
  binID=`awk -F '\t' -v orfFastaID="$orfFastaID" '$1 == orfFastaID { print $2 }' ORFs_G2B.tsv`
  python $HomeBio/AM3_binBasedAnalyses/doctor_petersons_neighborhood.py --bin_id $binID \
                                                                         --bin_file genomes/$binID.fna \
                                                                         --orf_fasta_id $orfFastaID \
                                                                         --ORF_location ORFs \
                                                                         --size_of_block 20000 \
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
echo "RUNNING KOFAMSCAN ON THE GENES IN THE NEIGHBORHOOD"
conda activate kofamscan
PYTHONPATH=""
PERL5LIB=""
analysis_name=prolix_hgcA_GN
input_orfs=$working_directory/GN_of_hgcA/hgcA_geneNeighborhood_all_orfs.faa
output_location=$working_directory/GN_of_hgcA/kofamscan_output
cd $working_directory
mkdir $output_location

rm -f $output_location/$analysis_name\_raw.tsv
echo "Running KOFAMscan for" $analysis_name
~/references/kofamscan_files/kofam_scan-1.3.0/exec_annotation -f detail-tsv \
                                                              -o $output_location/$analysis_name\_raw.tsv \
                                                              $input_orfs
conda deactivate



############################################
# Cluster ORFs with CD-HIT
############################################
echo "CLUSTERING THE GENES IN THE NEIGHBORHOOD"
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
# Pull all of these here:
# dataEdited/bins/binAnalysis/prolixibacteraceae_details/GN/GN_inputFiles

# 1. hgcA_geneNeighborhood_all_orfs_cluster_40.tsv
# 2. hgcA_id_list.txt
# 3. ORFs_G2B.tsv
# 4. prolix_hgcA_GN_raw.tsv
# Add data to spreadsheet: hgcA_geneNeighborhood_info.xlsx

# For KOFAM data, make sure to change it so that only the
# top annotation for each gene is there.

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
python $HomeBio/AM3_binBasedAnalyses/doctor_petersons_neighborhood_visualization.py --gff_file $working_directory/GN/hgcA_geneNeighborhood_all.gff \
                                                                                    --orf_data $working_directory/GN/hgcA_geneNeighborhood_info.txt \
                                                                                    --output_location $working_directory/GN/hgcA_geneNeighborhood_plot.pdf \
                                                                                    --pdf_height 10 \
                                                                                    --pdf_width 50

"""
