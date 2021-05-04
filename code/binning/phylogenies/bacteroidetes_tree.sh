#!/bin/sh

###########################
# code/binning/phylogenies/bacteroidetes_tree.sh
# Benjamin D. Peterson
###########################

cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood
mkdir phylogeny
mkdir phylogeny/Bacteroidetes
mkdir phylogeny/Bacteroidetes/ORFs
mkdir phylogeny/Bacteroidetes/lists


##################################################
##################################################
# Collect references
##################################################
##################################################

#########################
# Collect 5M Bacteroidetes hgcA+ bins
#########################
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/Bacteroidetes
grep -E 'BAC_00' ~/5M/dataEdited/binAnalysis/taxonomy/gtdbtk.bac120.summary.tsv | \
  awk -F '\t' '{ print $1 }' \
  > lists/5M_bacteroidetes_bin_list.txt
cat lists/5M_bacteroidetes_bin_list.txt | while read binID
do
  cp ~/5M/dataEdited/binAnalysis/ORFs/$binID.faa ORFs
done



#########################
# References selected for 5M project
#########################
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/Bacteroidetes
grep -E 'p__Bacteroidota' ~/5M/dataEdited/binAnalysis/phylogeny/bacteroidetes/reference_taxonomy.tsv | \
  grep -E 'RS_' | \
  awk -F '\t' '{ print $1"\t"$2 }' | \
  sed 's/RS_//' | \
  grep -v 'sp001027725' \
  > lists/reference_bacteroidetes_bin_list.txt
awk -F '\t' '{ print $1 }' lists/reference_bacteroidetes_bin_list.txt | while read binID
do
  cp -i ~/references/genomes/ORFs/$binID*.faa ORFs/$binID.faa
done


#########################
# hgcA+ references from McDaniel et al, 2020
#########################

# Done locally
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
    echo "Working on" $accessionID
    fileName=`ls $accessionID*`
    python $scripts/cleanFASTA.py $fileName
    mv $fileName\_temp.fasta $accessionID.fna
    rm $fileName
  else
    echo "Already cleaned" $accessionID
  fi
done

# Done on GLBRC
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/Bacteroidetes
mkdir scaffolds
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
scripts=~/HellsCanyon/code/generalUse
# Upload genomes to scaffolds folder, upload hgcA_bacteroidetes_genome_list.txt to lists folder
cat lists/hgcA_bacteroidetes_genome_list.txt | while read accessionID
do
  if [ ! -e  ORFs/$accessionID.faa ]; then
    prodigal -i scaffolds/$accessionID.fna \
              -o ORFs/$accessionID.gff \
              -f gff \
              -a ORFs/$accessionID.faa \
              -p single
    python $scripts/cleanFASTA.py ORFs/$accessionID.faa
    mv -f ORFs/$accessionID.faa_temp.fasta ORFs/$accessionID.faa
  else
    echo "Already predicted ORFs for" $accessionID
  fi
done



#########################
# Add bin from this project
#########################
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/
cp ORFs/anvio_hgcA_0130.faa phylogeny/Bacteroidetes/ORFs


#########################
# Get taxonomy for bins
#########################
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/Bacteroidetes
cat lists/reference_bacteroidetes_bin_list.txt


##################################################
##################################################
# Search for rp16 genes
##################################################
##################################################

#########################
# Concatenate all bins
#########################
scripts=~/HellsCanyon/code/generalUse/
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/Bacteroidetes/ORFs
cat *faa > ../ORFs.faa
$scripts/Fasta_to_Scaffolds2Bin.sh -e faa > ../ORFs_G2B.tsv

#########################
# Search for all genes
#########################
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/Bacteroidetes
mkdir rp16_hits

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
cat ~/references/rp16/rp16_list.txt | while read gene
do
  echo "Pulling" $gene "out"
  hmmsearch --tblout rp16_hits/$gene.out \
            --cpu 4 \
            --cut_tc \
            ~/references/rp16/$gene\_bact.HMM \
            ORFs.faa \
            > rp16_hits/$gene\_HMM_output.txt
  python $scripts/extract_protein_hitting_HMM.py \
                                    rp16_hits/$gene.out \
                                    ORFs.faa \
                                    rp16_hits/$gene.faa
  muscle -in rp16_hits/$gene.faa \
          -out rp16_hits/$gene.afa
done


#########################
# Rename fasta headers with bin name
#########################
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/Bacteroidetes
rm -f rp16_hits/G2B_key.tsv
grep -h '>' rp16_hits/*.afa | \
  sed 's/>//' | \
  while read geneID
do
  awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $0 }' ORFs_G2B.tsv >> rp16_hits/G2B_key.tsv
done
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/Bacteroidetes/rp16_hits
ls *.afa | while read binFile
do
  grep -h '>' $binFile | \
    sed 's/>//' | \
    while read geneID
  do
    binName=`awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $2 }' ../ORFs_G2B.tsv`
    sed -i "s/$geneID/$binName/" $binFile
  done
done



##################################################
##################################################
# Identify hgcA sequences in bins
##################################################
##################################################
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/Bacteroidetes
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
python $scripts/extract_protein_hitting_HMM.py \
                hgcA_hits/hgcA.out \
                ORFs.faa \
                hgcA_hits/hgcA_raw.faa
cd hgcA_hits
hmmalign -o hgcA_raw.sto \
            ~/references/hgcA/hgcA.hmm \
            hgcA_raw.faa
$scripts/convert_stockhold_to_fasta.py hgcA_raw.sto

# Make list of hgcA+ bins
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/Bacteroidetes
rm -f hgcA_bin_list.txt
grep '>' hgcA_hits/hgcA_raw.faa | sed 's/>//' | while read geneID
do
  awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $2 }' ORFs_G2B.tsv >> hgcA_hits/hgcA_bin_list.txt
done



##################################################
##################################################
# Generate tree
##################################################
##################################################

screen -S HCC_Bacteroidetes_tree
mkdir ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/Bacteroidetes/tree_building
# Upload alignment here
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/Bacteroidetes/tree_building
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
FastTree rp16_alignment_masked.afa > rp16.tree


screen -S HCC_Bacteroidetes_raxml_tree
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/Bacteroidetes/tree_building
raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s rp16_alignment_masked.afa \
        -n Bacteroidetes_rp16
