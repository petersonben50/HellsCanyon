#!/bin/sh

###########################
# code/binning/phylogenies/PVC_tree.sh
# Benjamin D. Peterson
###########################

cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood
mkdir phylogeny
mkdir phylogeny/PVC
mkdir phylogeny/PVC/ORFs
mkdir phylogeny/PVC/lists


##################################################
##################################################
# Collect references
##################################################
##################################################

#########################
# Collect 5M PVC hgcA+ bins
#########################
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/PVC
grep -E 'KIR_00|LEN_00' ~/5M/dataEdited/binAnalysis/taxonomy/gtdbtk.bac120.summary.tsv | \
  awk -F '\t' '{ print $1 }' \
  > lists/5M_PVC_bin_list.txt
cat lists/5M_PVC_bin_list.txt | while read binID
do
  cp ~/5M/dataEdited/binAnalysis/ORFs/$binID.faa ORFs
done


#########################
# Bins from Jones et al, 2019
#########################
scripts=~/HellsCanyon/code/generalUse/
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/PVC
grep 'PVC' ~/references/jonesGenomes/bin_names.tsv | \
  awk -F '\t' '{ print $1 }' | \
  while read binID
do
  ls ~/references/jonesGenomes/$binID.genes.faa
  python $scripts/cleanFASTA.py ~/references/jonesGenomes/$binID.genes.faa
  mv ~/references/jonesGenomes/$binID.genes.faa_temp.fasta ORFs/$binID.faa
done


#########################
# References selected for 5M project
#########################
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/PVC
grep -E 'c__Kiritimatiellae|c__Lentisphaeria' ~/5M/dataEdited/binAnalysis/phylogeny/PVC/reference_taxonomy.tsv | \
  awk -F '\t' '{ print $1"\t"$2 }' | \
  sed 's/GB_//' | \
  sed 's/RS_//' \
  > lists/reference_PVC_bin_list.txt
awk -F '\t' '{ print $1 }' lists/reference_PVC_bin_list.txt | while read binID
do
  cp ~/references/genomes/ORFs/$binID*.faa ORFs/$binID.faa
done


#########################
# Cultured isolates
#########################
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/PVC/ORFs
python $scripts/cleanFASTA.py GCF_003096415.1.faa
mv -f GCF_003096415.1.faa_temp.fasta GCF_003096415.1.faa

python $scripts/cleanFASTA.py GCF_900890425.1.faa
mv -f GCF_900890425.1.faa_temp.fasta GCF_900890425.1.faa

python $scripts/cleanFASTA.py GCF_900890705.1.faa
mv -f GCF_900890705.1.faa_temp.fasta GCF_900890705.1.faa


#########################
# Add bins from this project
#########################
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/
cp ORFs/anvio_hgcA_0220.faa phylogeny/PVC/ORFs
cp ORFs/anvio_hgcA_0261.faa phylogeny/PVC/ORFs
cp ORFs/anvio_hgcA_0040.faa phylogeny/PVC/ORFs
cp ORFs/anvio_hgcA_0110.faa phylogeny/PVC/ORFs


##################################################
##################################################
# Search for rp16 genes
##################################################
##################################################

#########################
# Concatenate all bins
#########################
scripts=~/HellsCanyon/code/generalUse/
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/PVC/ORFs
cat *faa > ../ORFs.faa
$scripts/Fasta_to_Scaffolds2Bin.sh -e faa > ../ORFs_G2B.tsv

#########################
# Search for all genes
#########################
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/PVC
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
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/PVC
rm -f rp16_hits/G2B_key.tsv
grep -h '>' rp16_hits/*.afa | \
  sed 's/>//' | \
  while read geneID
do
  awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $0 }' ORFs_G2B.tsv >> rp16_hits/G2B_key.tsv
done
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/PVC/rp16_hits
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
# Generate tree
##################################################
##################################################
screen -S HCC_PVC_tree
cd ~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/PVC/tree_building
raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s rp16_alignment_masked.afa \
        -n PVC_rp16
