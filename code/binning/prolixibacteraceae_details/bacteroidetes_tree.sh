#!/bin/sh

##################################################
##################################################

# code/binning/prolixibacteraceae_details/bacteroidetes_tree.sh
# Benjamin D. Peterson

# This script will generate the rp16-based tree for
# the Bacteroidetes genomes of interest from the
# study.
##################################################
##################################################
screen -S HCC_bacteroidetes
binAnalysisFolder=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis
HomeBio=/home/GLBRCORG/bpeterson26/HellsCanyon/code/HomeBio
cd $binAnalysisFolder
mkdir Bacteroidetes_tree
mkdir Bacteroidetes_tree/genomes
mkdir Bacteroidetes_tree/lists
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh



############################################
# Retrieve needed genomes
############################################

# Genomes from 5M project
cd $binAnalysisFolder/Bacteroidetes_tree
grep -E 'BAC_00' ~/5M/dataEdited/binAnalysis/taxonomy/gtdbtk.bac120.summary.tsv | \
  awk -F '\t' '{ print $1 }' \
  > lists/5M_bacteroidetes_bin_list.txt
cat lists/5M_bacteroidetes_bin_list.txt | while read binID
do
  cp ~/5M/dataEdited/binAnalysis/bins_processed/$binID.fna genomes
done

# References selected for 5M project
grep -E 'p__Bacteroidota' ~/5M/dataEdited/binAnalysis/phylogeny/bacteroidetes/reference_taxonomy.tsv | \
  grep -E 'RS_' | \
  awk -F '\t' '{ print $1 }' | \
  sed 's/RS_//' | \
  grep -v '001027725' \
  > lists/reference_bacteroidetes_bin_list.txt
cat lists/reference_bacteroidetes_bin_list.txt | while read binID
do
  if [ ! -e genomes/$binID.fna ]; then
    echo "Copy over" $binID
    cp ~/references/genomes/bins/$binID*.fna genomes/$binID.fna
  else
    echo $binID "already copied"
  fi
done

# hgcA+ references from McDaniel et al, 2020
# Done locally
"""
HCC
cd references/genomes/bacteroidetes
scripts=~/Documents/research/HellsCanyon/code/generalUse
ls *gz | while read file
do
  gunzip $file
done
# Clean fna and file name
list_of_interest=~/Documents/research/HellsCanyon/dataEdited/binning/phylogeny/bacteroidetes/hgcA_bacteroidetes_genome_list.txt
cat $list_of_interest | while read accessionID
do
  if [ ! -e $accessionID.fna ]; then
    echo 'Working on' $accessionID
    fileName=`ls $accessionID*`
    python $scripts/cleanFASTA.py $fileName
    mv $fileName\_temp.fasta $accessionID.fna
    rm $fileName
  else
    echo 'Already cleaned' $accessionID
  fi
done
# Upload all of these bins to /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/Bacteroidetes_tree/genomes
"""

# Add bins from this project
# Upload hgcA_minus_bins_to_use_for_phylogeny.txt to lists
cp ~/HellsCanyon/dataEdited/binning/bins_hgcA_keepers/DNA/anvio_hgcA_0130.fna genomes/anvio_hgcA_0130.fna
# Add hgcA- bins. Determined here: prolix_hgcA_minus_derep.R
cat lists/hgcA_minus_bins_to_use_for_phylogeny.txt | while read binID
do
  cp ~/HellsCanyon/dataEdited/binning/autoBinning/hqBinSet/DNA/$binID.fna genomes/$binID.fna
done

# Generate list of all bins to use for this tree.
ls genomes/*fna | \
    sed 's/genomes\///' | \
    sed 's/.fna//' > lists/bacteroidetes_bin_list.txt


############################################
# Predict ORFs
############################################
echo "PREDICTING ORFS WITH IMMA_ORF_STAN"
mkdir ORFs
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

# Predict ORFs
cat lists/bacteroidetes_bin_list.txt | while read accessionID
do
  if [ ! -e ORFs/$accessionID.faa ]; then
    echo "Predict ORFs for" $accessionID
    bash $HomeBio/PM4_binGeneration/IMMA_ORF_STAN_BINS.sh -b $accessionID \
                                                        -i genomes/$accessionID.fna \
                                                        -o ORFs \
                                                        -c $HomeBio/fasta_manipulation/cleanFASTA.py
  else
    echo "ORFs already predicted for" $accessionID
  fi
done

# Generate scaffold to bin file
rm -f ORFs_G2B.tsv ORFs.faa
bash $HomeBio/fasta_manipulation/Fasta_to_Scaffolds2Bin.sh -e faa \
                                                           -i ORFs \
                                                           > ORFs_G2B.tsv
cat ORFs/*.faa > ORFs.faa



############################################
# Pull out needed rp16 genes
############################################
cd $binAnalysisFolder/Bacteroidetes_tree
mkdir rp16_hits

cat ~/references/rp16/rp16_list.txt | while read gene
do
  echo "Pulling" $gene "out"
  hmmsearch --tblout rp16_hits/$gene.out \
            --cpu 4 \
            --cut_tc \
            ~/references/rp16/$gene\_bact.HMM \
            ORFs.faa \
            > rp16_hits/$gene\_HMM_output.txt
  python $HomeBio/fasta_manipulation/extract_protein_hitting_HMM.py \
                        rp16_hits/$gene.out \
                        ORFs.faa \
                        rp16_hits/$gene.faa
  muscle -in rp16_hits/$gene.faa \
          -out rp16_hits/$gene.afa
done


# Rename fasta headers with bin name
cd $binAnalysisFolder/Bacteroidetes_tree
rm -f rp16_hits/G2B_key.tsv
grep -h '>' rp16_hits/*.afa | \
  sed 's/>//' | \
  while read geneID
do
  awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $0 }' ORFs_G2B.tsv >> rp16_hits/G2B_key.tsv
done
cd $binAnalysisFolder/Bacteroidetes_tree/rp16_hits
ls *.afa | while read alignmentFile
do
  grep -h '>' $alignmentFile | \
    sed 's/>//' | \
    while read geneID
  do
    echo "Replacing" $geneID "with" $binName
    binName=`awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $2 }' ../ORFs_G2B.tsv`
    sed -i "s/$geneID/$binName/" $alignmentFile
  done
done



############################################
# Generate tree
############################################
cd $binAnalysisFolder/Bacteroidetes_tree
mkdir tree_building
# Upload alignment here
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
cd tree_building
FastTree rp16_alignment_masked.afa > rp16.tree


screen -S HCC_Bacteroidetes_raxml_tree
binAnalysisFolder=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis
cd $binAnalysisFolder/Bacteroidetes_tree/tree_building
raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s rp16_alignment_masked.afa \
        -n Bacteroidetes_rp16




############################################
# Look for hgcA
############################################
screen -S HCC_bacteroidetes
binAnalysisFolder=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis
HomeBio=/home/GLBRCORG/bpeterson26/HellsCanyon/code/HomeBio
cd $binAnalysisFolder/Bacteroidetes_tree
mkdir hgcA_hits
scripts=~/HellsCanyon/code/generalUse/
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

hmmsearch --tblout hgcA_hits/hgcA.out \
          --cpu 4 \
          --cut_tc \
          ~/references/hgcA/hgcA.hmm \
          ORFs.faa \
          > hgcA_hits/hgcA_HMM_output.txt
python $HomeBio/fasta_manipulation/extract_protein_hitting_HMM.py \
                hgcA_hits/hgcA.out \
                ORFs.faa \
                hgcA_hits/hgcA_raw.faa
muscle -in hgcA_hits/hgcA_raw.faa \
      -out hgcA_hits/hgcA_raw.afa


# Make list of hgcA+ bins
cd $binAnalysisFolder/Bacteroidetes_tree
rm -f hgcA_bin_list.txt
grep '>' hgcA_hits/hgcA_raw.faa | sed 's/>//' | while read geneID
do
  awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $2 }' ORFs_G2B.tsv >> hgcA_hits/hgcA_bin_list.txt
done



############################################
# Run references through GTDB database
############################################
conda deactivate
conda activate gtdbtk
PYTHONPATH=""
cd $binAnalysisFolder/Bacteroidetes_tree
rm -rf taxonomy
mkdir taxonomy

gtdbtk classify_wf \
        --cpus 16 \
        --extension fna \
        --genome_dir ./genomes \
        --out_dir taxonomy
awk -F '\t' '{ print $1"\t"$2 }' taxonomy/gtdbtk.*.summary.tsv \
        > taxonomy/taxonomy_summary.txt
conda deactivate
