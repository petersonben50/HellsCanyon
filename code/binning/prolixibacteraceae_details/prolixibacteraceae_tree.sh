#!/bin/sh


# Local computer
# Remove the GenBank entries for which there is a RefSeq entry
cd /Users/benjaminpeterson/Documents/research/HellsCanyon/dataEdited/bins/binAnalysis/phylogeny/prolixibacteraceae/reference_bins
ls GCF* | while read fileName
do
  GCF_ID=`echo $fileName | \
            sed 's/GCF_//' | \
            sed 's/_IMG.*//' | \
            sed 's/_ASM.*//'`
  rm GCA_$GCF_ID*
done

# On GLBRC
phylogeny=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny
scripts=~/HellsCanyon/code/generalUse
cd $phylogeny
mkdir prolixibacteraceae
mkdir prolixibacteraceae/reference_bins
# Upload zipped bins to reference_bins

# Unzip files
cd $phylogeny/prolixibacteraceae/reference_bins
ls *gz | while read file
do
  gunzip $file
done

# Clean up bins
ls *.fna | sed 's/\.1_.*/.1/' > accessionID_list.txt
cat accessionID_list.txt | while read accessionID
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


# Predicted ORFs
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
mkdir $phylogeny/prolixibacteraceae/ORFs
cd $phylogeny/prolixibacteraceae/ORFs
cat reference_bins/accessionID_list.txt | while read accessionID
do
  if [ ! -e  ORFs/$accessionID.faa ]; then
    prodigal -i reference_bins/$accessionID.fna \
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
# Add in outgroup bins and study bin
#########################
cd $phylogeny
cp Bacteroidetes/ORFs/anvio_hgcA_0130.faa prolixibacteraceae/ORFs/anvio_hgcA_0130.faa
cp Bacteroidetes/ORFs/GCF_000194605.1.faa prolixibacteraceae/ORFs/GCF_000194605.1.faa
cp Bacteroidetes/ORFs/GCF_000236705.1.faa prolixibacteraceae/ORFs/GCF_000236705.1.faa



##################################################
##################################################
# Add in hgcA- Prolixibacteraceae bins
##################################################
##################################################
screen -S HCC_Prolixibacteraceae
phylogeny=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny
hgcA_minus_bins=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binning/autoBinning/hqBinSet/ORFs
cd $phylogeny/prolixibacteraceae

cat hgcA_minus_bins_to_use_for_phylogeny.txt | while read accessionID
do
  echo "Copying over" $accessionID
  cp $hgcA_minus_bins/$accessionID.faa ORFs
done

##################################################
##################################################
# Search for rp16 genes
##################################################
##################################################

#########################
# Concatenate all bins
#########################
phylogeny=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny
scripts=~/HellsCanyon/code/generalUse/
cd $phylogeny/prolixibacteraceae
cat prolix_bins_to_remove.txt | while read accessionID
do
   rm -f ORFs/$accessionID*
done
cd ORFs
cat *faa > ../ORFs.faa
$scripts/Fasta_to_Scaffolds2Bin.sh -e faa > ../ORFs_G2B.tsv

#########################
# Search for all genes
#########################
cd $phylogeny/prolixibacteraceae
mkdir rp16_hits

screen -S HCC_prolix_rp16
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
phylogeny=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny
scripts=~/HellsCanyon/code/generalUse/
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
cd $phylogeny/prolixibacteraceae
rm -f rp16_hits/G2B_key.tsv
grep -h '>' rp16_hits/*.afa | \
  sed 's/>//' | \
  while read geneID
do
  awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $0 }' ORFs_G2B.tsv >> rp16_hits/G2B_key.tsv
done

cd $phylogeny/prolixibacteraceae/rp16_hits
ls *.afa | while read alignmentFile
do
  grep -h '>' $alignmentFile | \
    sed 's/>//' | \
    while read geneID
  do
    binName=`awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $2 }' ../ORFs_G2B.tsv`
    echo "Replacing" $geneID "with" $binName
    sed -i "s/$geneID/$binName/" $alignmentFile
  done
done



##################################################
##################################################
# Generate tree
##################################################
##################################################

screen -S HCC_Bacteroidetes_tree
phylogeny=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny
cd $phylogeny/prolixibacteraceae
mkdir tree_building
# Upload alignment here

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
cd tree_building
FastTree rp16_alignment_masked.afa > rp16.tree


screen -S HCC_Bacteroidetes_raxml_tree
phylogeny=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny
cd $phylogeny/prolixibacteraceae/tree_building
raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s rp16_alignment_masked.afa \
        -n prolixibacteraceae_rp16




##################################################
##################################################
# Run references through GTDB database
##################################################
##################################################


screen -S HCC_prolix_GTDB
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
conda activate gtdbtk
phylogeny=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny
cd $phylogeny/prolixibacteraceae

rm -rf taxonomy
mkdir taxonomy

gtdbtk classify_wf \
        --cpus 16 \
        --extension fna \
        --genome_dir ./reference_bins \
        --out_dir taxonomy
# Summarize them
cd taxonomy
awk -F '\t' '{ print $1"\t"$2 }' gtdbtk.*.summary.tsv \
        > taxonomy_summary.txt
conda deactivate




##################################################
##################################################
# Identify hgcA sequences in bins
##################################################
##################################################
cd $phylogeny/prolixibacteraceae
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
muscle -in hgcA_raw.faa \
      -out hgcA_raw.afa
# JABMQB010000125.1_1 has hgcA, but it's truncated right
# after the cap helix domain. We'll keep it for now.

# Make list of hgcA+ bins
cd $phylogeny/prolixibacteraceae
rm -f hgcA_bin_list.txt
grep '>' hgcA_hits/hgcA_raw.faa | sed 's/>//' | while read geneID
do
  awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $2 }' ORFs_G2B.tsv >> hgcA_hits/hgcA_bin_list.txt
done
