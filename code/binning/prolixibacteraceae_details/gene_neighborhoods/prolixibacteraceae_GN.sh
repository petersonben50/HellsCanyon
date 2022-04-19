#!/bin/sh





############################################
############################################
# Genomic context for hgcA
############################################
############################################

######################
# First pull out hgcA+ scaffolds
######################
screen -S HCC_prolix_hgcA
cd ~/HellsCanyon/dataEdited/binAnalysis
mkdir prolixibacteraceae_details
mkdir prolixibacteraceae_details/GN
cd prolixibacteraceae_details/GN
mkdir scaffolds
rm -f scaffolds/reference_scaffolds*

original_prolix_analysis=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny/prolixibacteraceae
grep '>' $original_prolix_analysis/hgcA_hits/hgcA_raw.faa | sed 's/>//' | while read geneID
do
  binID=`awk -F '\t' -v geneID="$geneID" '$1 == geneID { print $2 }' $original_prolix_analysis/ORFs_G2B.tsv`
  echo "working on bin" $binID", with gene" $geneID

  if [ $binID != "anvio_hgcA_0130" ]; then

    echo "Reference bin"
    scaffoldID=$(echo $geneID | awk -F '_' '{ print $1 }')
    grep -A 1 $scaffoldID$ $original_prolix_analysis/reference_bins/$binID.fna \
        >> scaffolds/reference_scaffolds.fna
    awk -v scaffoldID="$scaffoldID" '{ if ($1 == scaffoldID) print }' $original_prolix_analysis/ORFs/$binID.gff \
        >> scaffolds/reference_scaffolds.gff

  else
    echo "Study bin"
    scaffoldID=$(echo $geneID | cut -d '_' -f 1-2)
    assemblyName=$(echo $scaffoldID | awk -F '_' '{ print $1 }')
    grep -A 1 $scaffoldID$ ~/HellsCanyon/dataEdited/assemblies/scaffolds/$assemblyName\_assembly.fna \
        >> scaffolds/reference_scaffolds.fna

    awk -v scaffoldID="$scaffoldID" '{ if ($1 == scaffoldID) print }' ~/HellsCanyon/dataEdited/assemblies/ORFs/$assemblyName.gff \
        >> scaffolds/reference_scaffolds.gff

  fi
done



######################
# Isolate gene neighborhoods
######################

screen -S HCC_hgcA_gene_neighborhood
cd ~/HellsCanyon/dataEdited/binAnalysis/prolixibacteraceae_details/GN
mkdir GN_files
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate py_viz
PYTHONPATH=''
IFS=$'\n'
scripts=~/HellsCanyon/code/generalUse
original_prolix_analysis=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/phylogeny/prolixibacteraceae

grep '>' $original_prolix_analysis/hgcA_hits/hgcA_raw.faa | sed 's/>//' | while read hgcA_id
do
  if [ $hgcA_id == "HC18ME02_000000002532_3" ]; then

    scaffold_id=$(echo $hgcA_id | cut -d '_' -f 1-2)
    gene_id=$(echo $hgcA_id | \
                cut -d '_' -f 2-3 | \
                sed 's/^0*//g')

  else

    scaffold_id=$(echo $hgcA_id | cut -d '_' -f 1)
    ORF_ID=$(echo $hgcA_id | cut -d '_' -f 2)
    gene_id=$(awk -F '\t' -v scaffold_id="$scaffold_id" '$1 == scaffold_id { print $0 }' scaffolds/reference_scaffolds.gff | \
                  grep "ID=.*_$ORF_ID" | \
                  awk -F 'ID=' '{ print $2 }' | \
                  cut -d ';' -f 1)

  fi

  echo "Working on" $hgcA_id", from" $scaffold_id
  awk -F '\t' -v scaffold_id="$scaffold_id" '$1 == scaffold_id { print $0 }' scaffolds/reference_scaffolds.gff > GN_files/temp_scaffolds.gff
  echo "Searching for" $gene_id
  python $scripts/gene_neighborhood_extraction.py GN_files/temp_scaffolds.gff \
                                                  scaffolds/reference_scaffolds.fna \
                                                  $gene_id \
                                                  5000 \
                                                  GN_files/temp_$gene_id
  rm -f GN_files/temp_scaffolds.gff
done

cd GN_files
rm -f hgcA_geneNeighborhood.gff hgcA_geneNeighborhood.fna
cat temp_*.gff > hgcA_geneNeighborhood_raw.gff
cat temp_*.fna > hgcA_geneNeighborhood_raw.fna
rm -f *_neighborhood.*
