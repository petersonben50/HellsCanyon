#!/bin/sh

######################
# code/hgcA_analysis/hgcA_analysis.sh
# Benjamin D. Peterson

# This set of scripts will contain the bash
# code we need to pull out the hgcA sequences
# from our assembly and process them.
######################


############################################
############################################
# Identify hgcA sequences
############################################
############################################

######################
# Identify putative hgcA genes with HMM
######################

screen -S HCC_hgcA_search
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
scripts=~/HellsCanyon/code/generalUse/
ORFs=~/HellsCanyon/dataEdited/assemblies/ORFs
mkdir ~/HellsCanyon/dataEdited/hgcA_analysis
mkdir ~/HellsCanyon/dataEdited/hgcA_analysis/identification
cd ~/HellsCanyon/dataEdited/hgcA_analysis/identification
mkdir 2017 2018 2019
cd ~/HellsCanyon/dataEdited

cat ~/HellsCanyon/metadata/lists/assembly_key.csv | while read -r line
do
  assembly=`echo $line | awk -F ',' '{ print $1 }'`
  year=`echo $line | awk -F ',' '{ print $2 }'`
  echo "Searching for hgcA in" $assembly "from" $year
  if [ -e $ORFs/$assembly.faa ]; then
    if [ ! -e hgcA_analysis/identification/$year/$assembly\_hgcA_report.txt ]; then
      echo "Search for hgcA in" $assembly
      hmmsearch --tblout hgcA_analysis/identification/$year/$assembly\_hgcA.out \
                --cpu 4 \
                --cut_tc \
                ~/references/hgcA/hgcA.hmm \
                $ORFs/$assembly.faa \
                > hgcA_analysis/identification/$year/$assembly\_hgcA_report.txt
      python $scripts/extract_protein_hitting_HMM.py \
              hgcA_analysis/identification/$year/$assembly\_hgcA.out \
              $ORFs/$assembly.faa \
              hgcA_analysis/identification/$year/$assembly\_hgcA.faa
    else
      echo "Search is already done in" $assembly
    fi
  else
    echo "Genes aren't predicted for" $assembly
  fi
done


######################
# Concatenate and align all hgcA seqs for curation
######################
cd ~/HellsCanyon/dataEdited/hgcA_analysis/identification/
cat 201*/*_hgcA.faa > hgcA_raw.faa
# Align seqs
hmmalign -o hgcA_raw.sto \
            ~/references/hgcA/hgcA.hmm \
            hgcA_raw.faa
$scripts/convert_stockhold_to_fasta.py hgcA_raw.sto

# Curate the hgcA list, see notes in md file.
# Then:
cd ~/HellsCanyon/dataEdited/hgcA_analysis/identification/
grep '>' hgcA_good.afa | \
    sed 's/>//' \
    > hgcA_good.txt
sed 's/-//g' hgcA_good.afa > hgcA_good.faa


############################################
############################################
# Pull out depth of hgcA+ scaffolds
############################################
############################################

screen -S HCC_hgcA_depth
cd ~/HellsCanyon/dataEdited/hgcA_analysis
mkdir depth
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PERL5LIB=""
PYTHONPATH=""

awk -F ',' '{ print $1 }' ~/HellsCanyon/metadata/lists/metagenome_key.csv | while read metagenome
do
  cat identification/hgcA_good.txt | while read gene
  do
    scaffold=$(echo $gene | awk -F '_' '{ print $1"_"$2 }')
    assembly=$(echo $gene | awk -F '_' '{ print $1 }')
    if [ -e ~/HellsCanyon/dataEdited/mapping/$metagenome\_to_$assembly.bam ]; then
      echo "Calculating coverage of" $metagenome "over" $scaffold
      samtools depth -a -r $scaffold ~/HellsCanyon/dataEdited/mapping/$metagenome\_to_$assembly.bam \
          >> depth/$metagenome\_hgcA_depth_raw.tsv
    else
      echo $metagenome "not from same year as" $assembly "and" $gene "won't be found there"
    fi
  done

  echo "Aggregating hgcA depth information for" $metagenome
  python ~/HellsCanyon/code/generalUse/calculate_depth_length_contigs.py \
            depth/$metagenome\_hgcA_depth_raw.tsv \
            150 \
            depth/$metagenome\_hgcA_depth.tsv
  rm -f depth/$metagenome\_hgcA_depth_raw.tsv
done



############################################
############################################
# Genomic context for hgcA
############################################
############################################

######################
# First pull out hgcA+ scaffolds
######################
screen -S HCC_hgcA_context
cd ~/HellsCanyon/dataEdited/hgcA_analysis
# rm -r scaffolds
mkdir scaffolds

# Pull out scaffold FNA files
grep '>' identification/hgcA_raw.faa | \
    sed 's/>//' | \
    awk -F '_' '{ print $1"_"$2 }' \
    > scaffolds/hgcA_raw_scaffold_list.txt
rm scaffolds/hgcA_scaffolds.fna
cat scaffolds/hgcA_raw_scaffold_list.txt | while read scaffold
do
  assemblyName=$(echo $scaffold | awk -F '_' '{ print $1 }')
  echo "Pulling out" $scaffold "from" $assemblyName
  grep -A 1 $scaffold\$ ~/HellsCanyon/dataEdited/assemblies/scaffolds/$assemblyName\_assembly.fna \
      >> scaffolds/hgcA_scaffolds.fna
done

# Pull out GFF entries
rm -f scaffolds/hgcA_scaffolds.gff
cat scaffolds/hgcA_raw_scaffold_list.txt | while read scaffold
do
  assemblyName=$(echo $scaffold | awk -F '_' '{ print $1 }')
  echo "Pulling out" $scaffold "GFF entries"
  awk -v scaffold="$scaffold" '{ if ($1 == scaffold) print }' ~/HellsCanyon/dataEdited/assemblies/ORFs/$assemblyName.gff \
      >> scaffolds/hgcA_scaffolds.gff
done



######################
# Search downstream genes for hgcB
######################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate py_viz
PYTHONPATH=''
IFS=$'\n'
scripts=~/HellsCanyon/code/generalUse
ORFs=~/HellsCanyon/dataEdited/assemblies/ORFs
cd ~/HellsCanyon/dataEdited/hgcA_analysis
mkdir hgcB
rm -f hgcB/downstream_gene_list.txt
python $scripts/retrieve_downstream_gene_name.py \
          identification/hgcA_good.txt \
          scaffolds/hgcA_scaffolds.gff \
          hgcB/downstream_gene_list.txt

# Extract downstream amino acid sequences
rm -f hgcB/downstream_genes.faa
cat hgcB/downstream_gene_list.txt | while read gene
do
  assemblyName=$(echo $gene | cut -d "_" -f1)
  echo "Pulling out" $gene "faa entries"
  grep -A 1 $gene$ $ORFs/$assemblyName.faa >> hgcB/downstream_genes.faa
done
conda deactivate


# Search adjacent genes with hgcB HMM
cd ~/HellsCanyon/dataEdited/hgcA_analysis
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
scripts=~/HellsCanyon/code/generalUse
hmmsearch --tblout hgcB/hgcB.out \
          --cpu 4 \
          -T 30 \
          ~/references/hgcA/hgcB_5M.HMM \
          hgcB/downstream_genes.faa \
          > hgcB/hgcB_report.txt

# Check the hgcB hits
rm -f hgcB/hgcB.faa
grep -v "#" hgcB/hgcB.out | awk '{ print $1 }' | while read geneID
do
  grep -A 1 $geneID$ hgcB/downstream_genes.faa >> hgcB/hgcB.faa
done
grep '>' hgcB/hgcB.faa | wc -l
# Align sequences to HMM
hmmalign -o hgcB/hgcB.sto \
            ~/references/hgcA/hgcB_5M.HMM \
            hgcB/hgcB.faa
# Convert alignment to fasta format
$scripts/convert_stockhold_to_fasta.py hgcB/hgcB.sto


# Check sequences in Geneious.
# They all check out, so I'll just use the hgcB.faa
# file as our final set of hgcB seqs.
cd ~/HellsCanyon/dataEdited/hgcA_analysis/hgcB
grep '>' hgcB.faa | \
  sed 's/>//' > hgcB.txt

# Pull out non-hgcB downstream genes
cd ~/HellsCanyon/dataEdited/hgcA_analysis/hgcB
cat downstream_gene_list.txt hgcB.txt | \
  sort | \
  uniq -u
grep -A 1 'fall2017cluster6_000000000428_49' downstream_genes.faa
grep -A 1 'fall2017coassembly_000000448578_3' downstream_genes.faa
grep -A 1 'fall2017coassembly_000001334838_1' downstream_genes.faa
grep -A 1 'HC18HY300_000000013755_6' downstream_genes.faa
grep -A 1 'HC18ME02_000000018262_3' downstream_genes.faa
grep -A 1 'KMBP004F_000000216801_1' downstream_genes.faa
grep -A 1 'KMBP009B_000000084934_2' downstream_genes.faa




######################
# Isolate gene neighborhoods
######################

screen -S HCC_hgcA_gene_neighborhood
cd ~/HellsCanyon/dataEdited/hgcA_analysis
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate py_viz
PYTHONPATH=''
IFS=$'\n'
scripts=~/HellsCanyon/code/generalUse

grep '>' identification/hgcA_raw.faa | \
    sed 's/>//' | \
    while read hgcA_id
    do
      echo "Working on" $hgcA_id
      scaffold_id=$(echo $hgcA_id | cut -d '_' -f 1-2)
      awk -F '\t' -v scaffold_id="$scaffold_id" '$1 == scaffold_id { print $0 }' scaffolds/hgcA_scaffolds.gff > scaffolds/temp_scaffolds.gff
      gene_id=$(echo $hgcA_id | \
                  cut -d '_' -f 2-3 | \
                  sed 's/^0*//g')
      echo "Searching for" $gene_id
      python $scripts/gene_neighborhood_extraction.py scaffolds/temp_scaffolds.gff \
                                                      scaffolds/hgcA_scaffolds.fna \
                                                      $gene_id \
                                                      5000 \
                                                      scaffolds/temp_$gene_id

      rm -f scaffolds/temp_scaffolds.gff
    done

cd scaffolds
rm -f hgcA_geneNeighborhood.gff hgcA_geneNeighborhood.fna
cat temp_*.gff > hgcA_geneNeighborhood_raw.gff
cat temp_*.fna > hgcA_geneNeighborhood_raw.fna
rm -f *_neighborhood.*


############################################
############################################
# Classify hgcA seqs with pplacer workflow
############################################
############################################

screen -S HCC_hgcA_pplacer
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""
references=~/HellsCanyon/references/hgcA
workingDirectory=~/HellsCanyon/dataEdited/hgcA_analysis
scripts=~/HellsCanyon/code/generalUse
mkdir $workingDirectory/classification

# Generate fasta file of reference alignment
cp -avr ~/Everglades/references/hgcA/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg $references/
cd $references/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg
$scripts/convert_stockhold_to_fasta.py Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.stockholm
mv Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.afackholm $workingDirectory/classification/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.afa

# Generate alignment of sequences of interest
cd $workingDirectory
muscle -in identification/hgcA_good.faa \
        -out classification/hgcA_muscle.afa
cd classification/

# Combine the alignments of seqs from this study and references
muscle -profile -in1 hgcA_muscle.afa \
        -in2 Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.afa \
        -out hgcA_for_classification.afa
# Convert to stockholm format
python $scripts/convert_fasta_to_stockholm.py hgcA_for_classification.afa
conda deactivate

# Run pplacer
conda activate hgcA_classifier
pplacer -p \
        --keep-at-most 1 \
        --max-pend 1 \
        -c $references/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg/ \
        hgcA_for_classification.sto

# Make sqlite database of reference
rppr prep_db \
      --sqlite Hg_MATE_classify \
      -c $references/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg

# Generate taxonomic assignments using guppy
guppy classify -c $references/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg \
               --pp \
               --sqlite Hg_MATE_classify \
               hgcA_for_classification.jplace

# Save out this data to csv
guppy to_csv --point-mass \
              --pp \
              -o hgcA_classification.csv \
              hgcA_for_classification.jplace

# Visualize placements on FastTree
guppy tog --pp \
          -o hgcA_classification.nwk \
          hgcA_for_classification.jplace





############################################
############################################
# hgcA dereplication
############################################
############################################

cd ~/HellsCanyon/dataEdited/hgcA_analysis
mkdir dereplication
cdhit=~/programs/cdhit-master
$cdhit/cd-hit -g 1 \
              -i identification/hgcA_good.faa \
              -o dereplication/hgcA_good_acrossYear.faa \
              -c 0.97 \
              -n 5 \
              -d 0
$cdhit/clstr2txt.pl dereplication/hgcA_good_acrossYear.faa.clstr \
  > dereplication/hgcA_good_acrossYear.tsv



############################################
############################################
# Phylogenetic analysis of hgcA
############################################
############################################
screen -S HCC_hgcA_tree
cd ~/HellsCanyon/dataEdited/hgcA_analysis
workingDirectory=~/HellsCanyon/dataEdited/hgcA_analysis/phylogeny
scripts=~/HellsCanyon/code/generalUse
references=~/HellsCanyon/references/hgcA
mkdir $workingDirectory

# Pull out sequences for phylogenetic analysis.
# Add hgcA_rep_list.txt to the phylogeny folder
cd ~/HellsCanyon/dataEdited/hgcA_analysis/
rm -f phylogeny/hgcA_for_phylogeny.faa
cat phylogeny/hgcA_rep_list.txt | while read hgcA
do
  echo "Including" $hgcA "in phylogenetic analysis"
  grep -A 1 $hgcA identification/hgcA_good.faa >> phylogeny/hgcA_for_phylogeny.faa
done

# Generate alignment
cd $workingDirectory
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=''
muscle -in hgcA_for_phylogeny.faa \
        -out hgcA_for_phylogeny.afa
muscle -profile \
        -in1 hgcA_for_phylogeny.afa \
        -in2 $references/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full_align-bmge.fasta \
        -out hgcA_for_phylogeny_raw.afa

# Generate rough tree
FastTree hgcA_for_phylogeny_raw.afa \
    > rough_hgcA.tree
      # Download this to my local computer.

# Upload the list of sequences to remove
python $scripts/remove_fasta_seqs_using_list_of_headers.py hgcA_for_phylogeny_raw.afa \
                                                            seqs_to_remove.txt \
                                                            hgcA_for_phylogeny.afa

#########################
# Generate tree in RAxML
#########################

screen -S HCC_hgcA_tree
cd ~/HellsCanyon/dataEdited/hgcA_analysis/phylogeny
raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
#$raxml -f a \
#        -p 283976 \
#        -m PROTGAMMAAUTO \
#        -N autoMRE \
#        -x 2381 \
#        -T 20 \
#        -s hgcA_for_phylogeny_masked.afa \
#        -n hgcA
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s hgcA_for_phylogeny_masked.afa.reduced \
        -n hgcA



#########################
# Generate good tree with subset of references using RAxML
#########################

# Done locally
HCC
cp /Users/benjaminpeterson/Documents/research/Hg_MATE/versions/v1.01142021/Hg-MATE-Db.v1.01142021_ISOCELMAG_Hgc.fas \
    references
mkdir dataEdited/hgcA_analysis/phylogeny/final
mkdir dataEdited/hgcA_analysis/phylogeny/final/refs
finalFolder=dataEdited/hgcA_analysis/phylogeny/final/refs
rm -f $finalFolder/HgMate_reference_seqs_to_use.faa
cut -d"_" -f1-3 $finalFolder/reference_names_to_use.txt | while read reference_name
do
  grep -A 1 \
    $reference_name \
    references/Hg-MATE-Db.v1.01142021_ISOCELMAG_Hgc.fas \
    >> $finalFolder/HgMate_reference_seqs_to_use.faa
done
grep ">" $finalFolder/HgMate_reference_seqs_to_use.faa | \
  sed 's/>//' | tr -d '[:blank:]' \
  > $finalFolder/IDed_seqs.txt
wc -l $finalFolder/*
#cat $finalFolder/IDed_seqs.txt $finalFolder/reference_names_to_use.txt | \
#  sort

# Done on GLBRC
screen -S HCC_hgcA_tree
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=''
mkdir ~/HellsCanyon/dataEdited/hgcA_analysis/phylogeny/final
cd ~/HellsCanyon/dataEdited/hgcA_analysis/phylogeny/final
# Upload needed references
cat 5M_bin_seqs.faa \
    HgMate_reference_seqs_to_use.faa \
    jones_hgcA_seqs.faa \
    ../hgcA_for_phylogeny.faa \
    hgcA_paralogs_for_rooting.faa \
    > hgcA_for_tree_final.faa

# Generate alignment
muscle -in hgcA_for_tree_final.faa \
        -out hgcA_for_tree_final.afa

# Upload masked alignment (50% gaps)
# Then run RAxML to generate tree
raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s hgcA_for_tree_final_masked.afa \
        -n hgcA



#########################
# Generate tree in RAxML for UNIFRAC
#########################
screen -S UNIFRAC_hgcA
cd ~/HellsCanyon/dataEdited/hgcA_analysis/phylogeny
mkdir for_unifrac
cat hgcA_rep_list.txt | while read hgcA
do
  echo "Including" $hgcA "in phylogenetic analysis"
  grep -A 1 $hgcA hgcA_for_phylogeny.faa >> for_unifrac/hgcA.faa
done
cd for_unifrac
muscle -in hgcA.faa \
        -out hgcA.afa
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
trimal -in hgcA.afa \
       -out hgcA_trimmed.afa \
       -gt 0.5

raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s hgcA_trimmed.afa \
        -n hgcA
