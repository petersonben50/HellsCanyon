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
# Gene neighborhood analysis: Set up gene neighborhood analysis
##################################################
echo "PULLING OUT GENE NEIGHBORHOOD INFORMATION FOR HGCA GENES"

conda activate py_viz
PYTHONPATH=''
IFS=$'\n'
mkdir $working_directory/hgcA_analysis/GN_of_hgcA
cd $working_directory

cat hgcA_analysis/PVC_hgcA_geneID_list.txt | while read orfFastaID
do
  binID=`awk -F '\t' -v orfFastaID="$orfFastaID" '$1 == orfFastaID { print $2 }' HCC_PVC_ORFs_G2B.tsv`
  python $HomeBio/AM3_binBasedAnalyses/doctor_petersons_neighborhood.py --bin_id $binID \
                                                                         --bin_file genomes/$binID.fna \
                                                                         --orf_fasta_id $orfFastaID \
                                                                         --ORF_location ORFs \
                                                                         --size_of_block 20000 \
                                                                         --outputLocation hgcA_analysis/GN_of_hgcA
done

# Concatenate the output files
cd hgcA_analysis/GN_of_hgcA
rm -f hgcA_geneNeighborhood_all_orfs.faa
cat *_neighborhood.gff > hgcA_geneNeighborhood_all.gff
cat *_neighborhood.fna > hgcA_geneNeighborhood_all.fna
cat *.faa > hgcA_geneNeighborhood_all_orfs.faa
rm -f *_[0-9]*.faa
rm -f *_neighborhood.*
conda deactivate
grep '>' hgcA_geneNeighborhood_all_orfs.faa | \
  sed 's/>//' \
  > hgcA_geneNeighborhood_all_list.txt


############################################
# Gene neighborhood analysis: Run KOFAMSCAN on the neighborhood ORFs
############################################
echo "RUNNING KOFAMSCAN ON THE GENES IN THE NEIGHBORHOOD"
conda activate kofamscan
PYTHONPATH=""
PERL5LIB=""
cd $working_directory/hgcA_analysis/GN_of_hgcA
mkdir kofamscan_output

rm -f $working_directory/hgcA_analysis/GN_of_hgcA/kofamscan_output/$analysis_name\_raw.tsv
echo "Running KOFAMscan for" $analysis_name
~/references/kofamscan_files/kofam_scan-1.3.0/exec_annotation -f detail-tsv \
                                                              -o $working_directory/hgcA_analysis/GN_of_hgcA/kofamscan_output/GN_kofam_raw.tsv \
                                                              hgcA_geneNeighborhood_all_orfs.faa
conda deactivate



############################################
# Gene neighborhood analysis: Cluster ORFs with CD-HIT
############################################
echo "CLUSTERING THE GENES IN THE NEIGHBORHOOD"
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""

cd $working_directory/hgcA_analysis/GN_of_hgcA
mkdir clustering
cd-hit -g 1 \
        -i hgcA_geneNeighborhood_all_orfs.faa \
        -o clustering/hgcA_geneNeighborhood_all_orfs_cluster_40.faa \
        -c 0.40 \
        -n 2 \
        -d 0
clstr2txt.pl clustering/hgcA_geneNeighborhood_all_orfs_cluster_40.faa.clstr \
  > clustering/hgcA_geneNeighborhood_all_orfs_cluster_40.tsv



############################################
# Gene neighborhood analysis: Search for hgcB
############################################
cd $working_directory/hgcA_analysis/GN_of_hgcA
mkdir hgcB
python $HomeBio/fasta_manipulation/retrieve_downstream_genes.py \
        --gene_list_name $working_directory/hgcA_analysis/PVC_hgcA_geneID_list.txt \
        --gff_file $working_directory/HCC_PVC_ORFs.gff \
        --orf_file $working_directory/hgcA_analysis/GN_of_hgcA/hgcA_geneNeighborhood_all_orfs.faa \
        --output_file $working_directory/hgcA_analysis/GN_of_hgcA/hgcB/downstream_genes.faa
python $HomeBio/bin/protein_identification_and_alignment.py \
    --analysis_name PVC_hgcB \
    --hmm_file ~/references/hgcA/hgcB_5M.HMM \
    --orf_file $working_directory/hgcA_analysis/GN_of_hgcA/hgcB/downstream_genes.faa \
    --output_location $working_directory/hgcA_analysis/GN_of_hgcA/hgcB \
    --cpus_to_use 10 \
    --hmm_cutoff 30
grep '>' hgcB/PVC_hgcB.faa | \
    sed 's/>//' > hgcB/PVC_hgcB_list.txt


############################################
# Gene neighborhood analysis: Manually check the hgcB-sized seqs.
############################################
cd $working_directory/hgcA_analysis/GN_of_hgcA/hgcB
echo -e 'KIR_0001_19_36\nLEN_0037_131_5\nDHTZ01000076.1_16' | while read geneID
do
  grep -A 1 $geneID $working_directory/hgcA_analysis/GN_of_hgcA/hgcB/downstream_genes.faa >> hgcB_suspects.faa
done
cat hgcB_suspects.faa PVC_hgcB.faa > hgcB_suspects_interrogation.faa
muscle -align hgcB_suspects_interrogation.faa -output hgcB_suspects_interrogation.afa


############################################
# Figure for gene neighborhood
############################################
#
# Pull all of the listed files here:
# dataEdited/bins/binAnalysis/PVC_details/hgcA_analysis/GN_of_hgcA/GN_inputFiles

# 1. hgcA_geneNeighborhood_all_orfs_cluster_40.tsv
# 2. PVC_hgcA_geneID_list.txt
# 3. HCC_PVC_ORFs_G2B.tsv
# 4. GN_kofam_raw.tsv
# 5. hgcA_geneNeighborhood_all_list.txt

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

hgcA_tree=/Users/benjaminpeterson/Documents/research/HellsCanyon/results/bins/binAnalysis/PVC_details/hgcA_tree.pdf
working_directory=/Users/benjaminpeterson/Documents/research/HellsCanyon/dataEdited/bins/binAnalysis/PVC_details/hgcA_analysis/GN_of_hgcA
HomeBio=/Users/benjaminpeterson/Documents/programs/HomeBio
python $HomeBio/AM3_binBasedAnalyses/doctor_petersons_neighborhood_visualization.py --gff_file $working_directory/hgcA_geneNeighborhood_all.gff \
                                                                                    --orf_data $working_directory/hgcA_geneNeighborhood_info_PVC.txt \
                                                                                    --output_location $working_directory/hgcA_geneNeighborhood_plot.pdf \
                                                                                    --gene_tree_file $hgcA_tree \
                                                                                    --gn_x_start 15000 \
                                                                                    --gn_x_end 26000


"""
