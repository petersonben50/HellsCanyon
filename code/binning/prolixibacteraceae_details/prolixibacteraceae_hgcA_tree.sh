#!/bin/sh

conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

working_directory=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/prolixibacteraceae_details
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh

cd $working_directory
mkdir hgcA_tree

##########################
# Bring in hgcA genes
##########################
original_phylogeny_analysis=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny
cp $original_phylogeny_analysis/prolixibacteraceae/hgcA_hits/hgcA_raw.faa hgcA_tree/hgcA_prolix.faa
# Dereplicate them. Basically removes the study sequences
cd-hit-2d -i hgcA_tree/hgcA_prolix.faa \
          -i2 $original_phylogeny_analysis/Bacteroidetes/hgcA_hits/hgcA_raw.faa \
          -o hgcA_tree/hgcA_bacteroidetes.faa \
          -c 0.97 \
          -n 5
cat hgcA_tree/hgcA_prolix.faa hgcA_tree/hgcA_bacteroidetes.faa > hgcA_tree/hgcA_bact.faa
rm -f hgcA_tree/hgcA_prolix.faa hgcA_tree/hgcA_bacteroidetes.faa


##########################
# Get bin data
##########################
cd hgcA_tree
grep '>' hgcA_bact.faa | sed 's/>//' | while read geneID
do
  awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $0 }' ~/HellsCanyon/dataEdited/binAnalysis/phylogeny/Bacteroidetes/ORFs_G2B.tsv >> hgcA_G2B.tsv
done

grep '>' hgcA_bact.faa | sed 's/>//' | while read geneID
do
  if [ $geneID != "HC18ME02_000000002532_3" ];then
    awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $0 }' ~/HellsCanyon/dataEdited/binAnalysis/prolixibacteraceae_details/ORFs_G2B.tsv >> hgcA_G2B.tsv
  fi
done



##########################
# Prepare alignment for tree
##########################
muscle -in hgcA_bact.faa \
        -out hgcA_bact.afa
        # Trim the alignment
trimal -in hgcA_bact.afa \
        -out hgcA_bact_trimmed.afa \
        -gt 0.5



##########################
# Generate tree
##########################
FastTree hgcA_bact_trimmed.afa > hgcA_bact.tree
raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s hgcA_bact_trimmed.afa \
        -n hgcA_bact
