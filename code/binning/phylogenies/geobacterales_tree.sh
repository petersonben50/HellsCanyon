#!/bin/sh

###########################
# code/binning/phylogenies/geobacterales_tree.sh
# Benjamin D. Peterson
###########################


######################################################
# Get bins using E-Utilities and wget
######################################################
# On local computer
# Separate into RefSeq and GenBank sequences
cd /Users/benjaminpeterson/Documents/research/HellsCanyon/dataEdited/bins/binAnalysis/phylogeny/geobacterales
grep "GCF_" GTDB_ref_accession_numbers.txt > GTDB_RefSeq_accession_numbers.txt
grep "GCA_" GTDB_ref_accession_numbers.txt > GTDB_GenBank_accession_numbers.txt

# E-Utilities to get the FTP links.
epost -db assembly \
      -input GTDB_RefSeq_accession_numbers.txt \
      -format acc | \
  esummary | \
  xtract -pattern DocumentSummary -element FtpPath_RefSeq \
  > ftp_paths.txt
epost -db assembly \
      -input GTDB_GenBank_accession_numbers.txt \
      -format acc | \
  esummary | \
  xtract -pattern DocumentSummary -element FtpPath_GenBank \
  >> ftp_paths.txt

# Retrieve files with wget
cd NCBI_downloads
wget -r -nd -i ../ftp_paths.txt

# Pull out the DNA files
cd ..
mkdir DNA
cp NCBI_downloads/*genomic.fna* DNA
cd DNA
rm -f *_cds_from_genomic* *_rna_from_genomic*


######################################################
# Clean up genomic bins
######################################################
screen -S HCC_geobacterales
phylogeny=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny
scripts=~/HellsCanyon/code/generalUse
cd $phylogeny
mkdir geobacterales
mkdir geobacterales/reference_bins
# Upload the reference genome files to GLBRC
# Also upload GTDB_ref_accession_numbers.txt
cd geobacterales/reference_bins
ls *gz | while read file
do
  gunzip $file
done
grep -v "anvio" ../GTDB_ref_accession_numbers.txt | while read accessionID
do
  if [ ! -e $accessionID.fna ]; then
    echo "Working on" $accessionID
    fileName=`ls $accessionID*`
    python $scripts/cleanFASTA.py $fileName
    mv $fileName\_temp.fasta $accessionID.fna
    #rm $fileName
  else
    echo "Already cleaned" $accessionID
  fi
done


######################################################
# Identify proteins from bins
######################################################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
phylogeny=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny
scripts=~/HellsCanyon/code/generalUse
cd $phylogeny/geobacterales
mkdir ORFs

grep -v "anvio" GTDB_ref_accession_numbers.txt | while read accessionID
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

# Add in study bin
cd $phylogeny
cp ~/HellsCanyon/dataEdited/binning/bins_hgcA/ORFs/anvio_hgcA_0210.faa geobacterales/ORFs/
# Add in bin from Mendota
cp ~/5M/dataEdited/binAnalysis/ORFs/GEO_0030.faa geobacterales/ORFs/

#########################
# Concatenate all bins
#########################
phylogeny=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny
scripts=~/HellsCanyon/code/generalUse/
cd $phylogeny/geobacterales/ORFs
cat *faa > ../ORFs.faa
$scripts/Fasta_to_Scaffolds2Bin.sh -e faa > ../ORFs_G2B.tsv



##################################################
##################################################
# Search for rp16 genes
##################################################
##################################################

#########################
# Search for all genes
#########################

screen -S HCC_geobac_rp16
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
phylogeny=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny
cd $phylogeny/geobacterales
mkdir rp16_hits
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
cd $phylogeny/geobacterales
rm -f rp16_hits/G2B_key.tsv
grep -h '>' rp16_hits/*.afa | \
  sed 's/>//' | \
  while read geneID
do
  awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $0 }' ORFs_G2B.tsv >> rp16_hits/G2B_key.tsv
done

cd $phylogeny/geobacterales/rp16_hits
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
cd $phylogeny/geobacterales
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
cd $phylogeny/geobacterales/tree_building
raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s rp16_alignment_masked.afa \
        -n geobacterales_rp16




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
cd $phylogeny/geobacterales

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
screen -S HCC_prolix_hgcA
phylogeny=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny
cd $phylogeny/geobacterales
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
cd $phylogeny/geobacterales
rm -f hgcA_bin_list.txt
grep '>' hgcA_hits/hgcA_raw.faa | sed 's/>//' | while read geneID
do
  awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $2 }' ORFs_G2B.tsv >> hgcA_hits/hgcA_bin_list.txt
done
