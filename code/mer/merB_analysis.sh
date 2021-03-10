#!/bin/sh

######################
# code/mer/merB_analysis.sh
# Benjamin D. Peterson

# This set of scripts will contain the bash
# code we need to pull out the merB sequences
# from our assembly and process them.
######################


############################################
############################################
# Identify merB sequences
############################################
############################################

######################
# Identify putative merB genes with HMM
######################

screen -S HCC_merB_search
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
scripts=~/HellsCanyon/code/generalUse/
ORFs=~/HellsCanyon/dataEdited/assemblies/ORFs
mkdir ~/HellsCanyon/dataEdited/merB_analysis
mkdir ~/HellsCanyon/dataEdited/merB_analysis/identification
cd ~/HellsCanyon/dataEdited/merB_analysis/identification
mkdir 2017 2018 2019
cd ~/HellsCanyon/dataEdited

cat ~/HellsCanyon/metadata/lists/assembly_key.csv | while read -r line
do
  assembly=`echo $line | awk -F ',' '{ print $1 }'`
  year=`echo $line | awk -F ',' '{ print $2 }'`
  if [ -e $ORFs/$assembly.faa ]; then
    if [ ! -e merB_analysis/identification/$year/$assembly\_merB_report.txt ]; then
      echo "Search for merB in" $assembly "from" $year
      hmmsearch --tblout merB_analysis/identification/$year/$assembly\_merB.out \
                --cpu 4 \
                --cut_tc \
                ~/references/merB/MerB.hmm \
                $ORFs/$assembly.faa \
                > merB_analysis/identification/$year/$assembly\_merB_report.txt
      python $scripts/extract_protein_hitting_HMM.py \
              merB_analysis/identification/$year/$assembly\_merB.out \
              $ORFs/$assembly.faa \
              merB_analysis/identification/$year/$assembly\_merB.faa
    else
      echo "Search is already done in" $assembly
    fi
  else
    echo "Genes aren't predicted for" $assembly
  fi
done

cd merB_analysis/identification
grep -v '#' */*_merB.out

# Get list of all merB proteins
grep '>' 201*/*_merB.faa | \
  sed 's/>//' | \
  cut -d":" -f2 \
  > merB_raw_list.txt


######################
# Concatenate and align all merB seqs for curation
######################
cd ~/HellsCanyon/dataEdited/merB_analysis/identification/
cat 201*/*_merB.faa > merB_raw.faa
# Align seqs
hmmalign -o merB_raw.sto \
            ~/references/merB/MerB.hmm \
            merB_raw.faa
$scripts/convert_stockhold_to_fasta.py merB_raw.sto


######################
# Concatenate and align all merB seqs with references for tree generation
######################
cat 201*/*_merB.faa \
    ../uniprot_ec_4.99.1.2.fasta \
    > merB_tree_all.faa
hmmalign -o merB_tree_all.sto \
            ~/references/merB/MerB.hmm \
            merB_tree_all.faa
$scripts/convert_stockhold_to_fasta.py merB_tree_all.sto
# Download alignment to computer, mask in Geneious
# at 50% gaps.

FastTree merB_tree_all_masked.afa \
    > merB_tree_for_ID.tree


######################
# Inspect individual clusters
######################
cd ~/HellsCanyon/dataEdited/merB_analysis/identification/
mkdir subclusters
# Upload lists.

python $scripts/cleanFASTA.py merB_tree_all.afa
mv -f merB_tree_all.afa_temp.fasta merB_tree_all.afa

ls subclusters/*_list.txt | while read file
do
  clusterID=`echo $file | sed 's/subclusters\///' | awk -F '_' '{ print $1 }'`
  rm -f subclusters/$clusterID.afa
  cat $file | while read seqID
  do
    grep -A 1 $seqID merB_tree_all.afa >> subclusters/$clusterID.afa
  done
done

rm -f subclusters/all_seqs_of_interest.afa
ls subclusters/*_list.txt | while read file
do
  clusterID=`echo $file | sed 's/subclusters\///' | awk -F '_' '{ print $1 }'`
  cat $file | while read seqID
  do
    grep -A 1 $seqID merB_tree_all.afa >> subclusters/all_seqs_of_interest.afa
  done
done
cd subclusters
FastTree all_seqs_of_interest_masked.afa > all_seqs_of_interest.tree


######################
# Make key for final groups
######################

cd ~/HellsCanyon/dataEdited/merB_analysis/identification/final_groups
echo "seqID,clusterName" > group_key.csv
ls *afa | while read file
do
  clusterName=`echo $file | sed 's/.afa//'`
  grep '>' $file | \
    sed 's/>//' | \
    awk -v clusterName="$clusterName" '{ print $0","clusterName}' \
    >> group_key.csv
done
grep 'KMBP\|fall2017\|HC18' group_key.csv > group_key_samplesOnly.csv

############################################
############################################
# Dereplicate sequences
############################################
############################################

cd ~/HellsCanyon/dataEdited/merB_analysis
mkdir dereplication
cdhit=~/programs/cdhit-master
awk -F ',' '{ print $2 }' ~/HellsCanyon/metadata/lists/assembly_key.csv | \
  sort | uniq | \
  while read year
  do
    cat identification/$year/*.faa > identification/$year/merB_from_$year.faa
    $cdhit/cd-hit -g 1 \
                  -i identification/$year/merB_from_$year.faa \
                  -o dereplication/merB_from_$year.faa \
                  -c 0.97 \
                  -n 5 \
                  -d 0
    $cdhit/clstr2txt.pl dereplication/merB_from_$year.faa.clstr \
      > dereplication/merB_from_$year.tsv
  done
cd dereplication
cat merB_from*.faa > merB_derep.faa
grep '>' merB_derep.faa | \
  sed 's/>//' \
  > merB_derep_list.txt

############################################
############################################
# Pull out depth of merB+ scaffolds
############################################
############################################

screen -S HCC_merB_depth
cd ~/HellsCanyon/dataEdited/merB_analysis
mkdir depth
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PERL5LIB=""
PYTHONPATH=""

awk -F ',' '{ print $1 }' ~/HellsCanyon/metadata/lists/metagenome_key.csv | while read metagenome
do
  cat identification/merB_raw_list.txt | while read gene
  do
    scaffold=$(echo $gene | awk -F '_' '{ print $1"_"$2 }')
    assembly=$(echo $gene | awk -F '_' '{ print $1 }')
    if [ -e ~/HellsCanyon/dataEdited/mapping/$metagenome\_to_$assembly.bam ]; then
      echo "Calculating coverage of" $metagenome "over" $scaffold
      samtools depth -a -r $scaffold ~/HellsCanyon/dataEdited/mapping/$metagenome\_to_$assembly.bam \
          >> depth/$metagenome\_merB_depth_raw.tsv
    else
      echo $metagenome "not from same year as" $assembly "and" $gene "won't be found there"
    fi
  done

  echo "Aggregating merB depth information for" $metagenome
  python ~/HellsCanyon/code/generalUse/calculate_depth_length_contigs.py \
            depth/$metagenome\_merB_depth_raw.tsv \
            150 \
            depth/$metagenome\_merB_depth.tsv
  rm -f depth/$metagenome\_merB_depth_raw.tsv
done




############################################
############################################
# Genomic context for merB
############################################
############################################

######################
# First pull out merB+ scaffolds
######################
screen -S HCC_merB_scaffolds
cd ~/HellsCanyon/dataEdited/merB_analysis
# rm -r scaffolds
mkdir scaffolds

# Identify scaffolds with multiple merB hits
awk -F '_' '{ print $1"_"$2 }' identification/merB_raw_list.txt | \
  sort | uniq -c | \
  awk '$1 == 2 { print $2 }' > scaffolds/scaffolds_with_multiple_merB_genes.txt
# Remove these from the neighborhood analysis,
# we'll look at them later.

# Pull out scaffold FNA files
awk -F '_' '{ print $1"_"$2 }' identification/merB_raw_list.txt | \
  sort | uniq -c | \
  awk '$1 == 1 { print $2 }' | \
  awk -F '_' '{ print $1"_"$2 }' \
  > scaffolds/scaffolds_for_neighborhood_analysis.txt

rm -f scaffolds/merB_scaffolds.fna
cat scaffolds/scaffolds_for_neighborhood_analysis.txt | while read scaffold
do
  assemblyName=$(echo $scaffold | awk -F '_' '{ print $1 }')
  echo "Pulling out" $scaffold "from" $assemblyName
  grep -A 1 $scaffold\$ ~/HellsCanyon/dataEdited/assemblies/scaffolds/$assemblyName\_assembly.fna \
      >> scaffolds/merB_scaffolds.fna
done

# Pull out GFF entries
rm -f scaffolds/merB_scaffolds.gff
cat scaffolds/scaffolds_for_neighborhood_analysis.txt | while read scaffold
do
  assemblyName=$(echo $scaffold | awk -F '_' '{ print $1 }')
  echo "Pulling out" $scaffold "GFF entries"
  awk -v scaffold="$scaffold" '{ if ($1 == scaffold) print }' ~/HellsCanyon/dataEdited/assemblies/ORFs/$assemblyName.gff \
      >> scaffolds/merB_scaffolds.gff
done



######################
# Isolate gene neighborhoods
######################
screen -S HCC_merB_neighborhoods
cd ~/HellsCanyon/dataEdited/merB_analysis
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate py_viz
PYTHONPATH=''
IFS=$'\n'
scripts=~/HellsCanyon/code/generalUse

# Genes for neighborhood analysis
grep -f scaffolds/scaffolds_for_neighborhood_analysis.txt \
        identification/merB_raw_list.txt \
        > scaffolds/genes_for_neighborhood_analysis.txt

cat scaffolds/genes_for_neighborhood_analysis.txt | while read merB_id
do
  echo "Working on" $merB_id
  scaffold_id=$(echo $merB_id | cut -d '_' -f 1-2)
  awk -F '\t' -v scaffold_id="$scaffold_id" '$1 == scaffold_id { print $0 }' scaffolds/merB_scaffolds.gff > scaffolds/temp_scaffolds.gff
  gene_id=$(echo $merB_id | \
              cut -d '_' -f 2-3 | \
              sed 's/^0*//g')
  echo "Searching for" $gene_id
  python $scripts/gene_neighborhood_extraction.py scaffolds/temp_scaffolds.gff \
                                                  scaffolds/merB_scaffolds.fna \
                                                  $gene_id \
                                                  10000 \
                                                  scaffolds/temp_$gene_id

  rm -f scaffolds/temp_scaffolds.gff
done

cd scaffolds
rm -f merB_geneNeighborhood.gff merB_geneNeighborhood.fna
cat temp_*.gff > merB_geneNeighborhood_raw.gff
cat temp_*.fna > merB_geneNeighborhood_raw.fna
rm -f *_neighborhood.*

######################
# Pull out amino acid sequences
######################
rm -f friendly_neighborhood_genes_list.txt
cat merB_geneNeighborhood_raw.gff | while read line
do
  scaffoldID=`echo $line | awk -F '\t' '{ print $1 }'`
  geneID=`echo $line | awk -F '\t' '{ print $9 }' | \
            sed 's/ID=//' | \
            cut -d"_" -f2 | \
            cut -d";" -f1`
  echo $scaffoldID"_"$geneID >> friendly_neighborhood_genes_list.txt
done
cd ~/HellsCanyon/dataEdited
cat merB_analysis/scaffolds/friendly_neighborhood_genes_list.txt | while read geneID
do
  assembly=`echo $geneID | awk -F '_' '{ print $1 }'`
  grep -A 1 $geneID$ assemblies/ORFs/$assembly.faa >> merB_analysis/scaffolds/friendly_neighborhood_genes.faa
done

#########################
# Run Kofamscan on ORFs
#########################
cp -avr ~/HeCaMT/code/kofam_scan-1.3.0 ~/HellsCanyon/code/
# config.yml file updated as part of HeCaMT analysis

screen -S kofamscan
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate kofamscan
PYTHONPATH=""
PERL5LIB=""
cd ~/HellsCanyon/code/kofam_scan-1.3.0
scaffolds=~/HellsCanyon/dataEdited/merB_analysis/scaffolds
mkdir $scaffolds/KOFAM_output

./exec_annotation -f detail-tsv \
                  -o $scaffolds/KOFAM_output/friendly_neighborhood_genes_annotationsAll.tsv \
                  $scaffolds/friendly_neighborhood_genes.faa
cd $scaffolds/KOFAM_output
grep '*' friendly_neighborhood_genes_annotationsAll.tsv > friendly_neighborhood_genes_annotations.tsv





#########################
# Generate gene neighborhood images
#########################

HCC
cd dataEdited/mer/scaffolds
mkdir neighborhood_viz
scripts=~/Documents/research/HellsCanyon/code/mer
keyFile=~/Documents/research/HellsCanyon/dataEdited/mer/identification/final_groups/group_key_samplesOnly.csv

cat $keyFile | while read line
do
  gene_of_interest=$(echo $line | awk -F ',' '{ print $1 }')
  scaffold=$(echo $gene_of_interest | rev | cut -d"_" -f2- | rev)
  group_of_interest=$(echo $line | awk -F ',' '{ print $2 }')

  echo "Working on" $gene_of_interest "in group" $group_of_interest
  awk -F '\t' -v "scaffold=$scaffold" '$1 == scaffold { print $0 }' merB_geneNeighborhood_annot.gff \
    > temp.gff
  python $scripts/extract_viz_table_from_gff_merB.py \
      temp.gff \
      $gene_of_interest \
      neighborhood_viz/$gene_of_interest\_from_$group_of_interest\_scaffold.csv
  rm temp.gff
done

cd neighborhood_viz/
awk -F ',' '{ print $2 }' $keyFile | \
  sort | uniq | \
  while read group_of_interest
  do
    echo "Working on" $group_of_interest
    fileForHeader=$(ls *_from_$group_of_interest\_scaffold.csv | head -n 1)
    head -n 1 $fileForHeader > $group_of_interest\_forFigure.csv
    for file in *_from_$group_of_interest\_scaffold.csv
    do
      tail -n +2 $file >> $group_of_interest\_forFigure.csv
    done
  done

rm -f *_scaffold.csv


# See protocol for details
conda activate py_viz
#results=/Users/benjaminpeterson/Documents/research/HellsCanyon/results/
awk -F ',' '{ print $2 }' $keyFile | \
  sort | uniq | \
  while read group_of_interest
  do
    python $scripts/merB_geneNeighborhood_figure.py \
                $group_of_interest\_forFigure.csv \
                $group_of_interest\_neighborhood.pdf
  done


conda deactivate
