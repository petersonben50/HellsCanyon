#!/bin/sh

#########################
# code/metabolic_analyses/metabolic_hmms.sh
# Benjamin D. Peterson
#########################


#########################
# Identification of metabolic proteins
#########################

screen -S HCC_metabolic_HMMs

cd ~/HellsCanyon/dataEdited/
mkdir metabolic_analyses
mkdir metabolic_analyses/identification

scripts=~/HellsCanyon/code/generalUse
metabolic_HMMs=~/HellsCanyon/references/metabolic_HMMs
ORFs=~/HellsCanyon/dataEdited/assemblies/ORFs

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

tail -n +2 $metabolic_HMMs.csv | awk -F ',' '{ print $2 }' | while read HMM
do
  cd ~/HellsCanyon/dataEdited/metabolic_analyses
  geneName=$(awk -F ',' -v HMM="$HMM" '$2 == HMM { print $1 }' $metabolic_HMMs.csv)

  if [ ! -e identification/$geneName.afa ]; then
    echo "Searching for" $geneName
    mkdir identification/$geneName

    cat ~/HellsCanyon/metadata/lists/assembly_key.csv | while read -r line
    do
      assembly=`echo $line | awk -F ',' '{ print $1 }'`
      year=`echo $line | awk -F ',' '{ print $2 }'`
      if [ ! -d identification/$geneName/$year ]; then
        mkdir identification/$geneName/$year
      fi

      if [ ! -e identification/$geneName/$year/$assembly.out ]; then
        echo "Searching for" $geneName "in" $assembly
        hmmsearch --tblout identification/$geneName/$year/$assembly\_$geneName.out \
                  --cpu 4 \
                  --cut_tc \
                  $metabolic_HMMs/$HMM \
                  $ORFs/$assembly.faa \
                  > identification/$geneName/$year/$assembly\_$geneName.txt
        lineCount=`wc -l < identification/$geneName/$year/$assembly\_$geneName.out`
        if [ $lineCount -eq 13 ]; then
          echo "No" $gene "hits in" $geneName
        else
          echo "Pulling" $geneName "sequences out of" $assembly
          python $scripts/extract_protein_hitting_HMM.py \
                  identification/$geneName/$year/$assembly\_$geneName.out \
                  $ORFs/$assembly.faa \
                  identification/$geneName/$year/$assembly\_$geneName.faa
        fi
      else
        echo "Search for" $geneName "is already done in" $assembly
      fi
    done

    # Aggregate all sequences and align to HMM
    cd identification
    awk -F ',' '{ print $2 }' ~/HellsCanyon/metadata/lists/assembly_key.csv | \
      sort | uniq | \
      while read year
      do
        shopt -s nullglob
        for i in $geneName/$year/*$geneName.faa; do FOUND=$i;break;done

        if [ ! -z $FOUND ]; then
          echo "Concatenating and aligning" $geneName
          cat $geneName/$year/*$geneName.faa \
              > $geneName\_from_$year\_all.faa
          hmmalign -o $geneName\_from_$year.sto \
                      $metabolic_HMMs/$HMM \
                      $geneName\_from_$year\_all.faa
          $scripts/convert_stockhold_to_fasta.py $geneName\_from_$year.sto
          grep '>' $geneName\_from_$year\_all.faa | \
            sed 's/>//' \
            > $geneName\_from_$year\_all_list.txt
        else
          echo "No" $geneName "sequences found at all :("
        fi
    done
    FOUND=""
  else
    echo "Already pulled out" $geneName "sequences"
  fi
done

conda deactivate

exit




#########################
# Extract depths of all scaffolds
#########################

# Generate list of all scaffolds that have potential metabolic genes on them
cd ~/HellsCanyon/dataEdited/metabolic_analyses
mkdir depth
grep '>' -h identification/201*/*.afa | \
  sed 's/>//' | \
  cut -d"_" -f1,2 | \
  sort | \
  uniq \
  > depth/scaffold_all_list.txt

# Pull out all depths
chmod +x /home/GLBRCORG/bpeterson26/HellsCanyon/code/executables/aggregate_depth_proteins.sh
condor_submit /home/GLBRCORG/bpeterson26/HellsCanyon/code/submission/aggregate_depth_proteins.sub


############################################
############################################
# Dereplicate sequences
############################################
############################################

cd ~/HellsCanyon/dataEdited/metabolic_analyses
mkdir dereplication
cdhit=~/programs/cdhit-master

echo "geneID,geneName" > metabolic_gene_key.csv

tail -n +2 $metabolic_HMMs.csv | awk -F ',' '{ print $1 }' | while read geneName
do
  cd ~/HellsCanyon/dataEdited/metabolic_analyses
  awk -F ',' '{ print $2 }' ~/HellsCanyon/metadata/lists/assembly_key.csv | \
    sort | uniq | \
    while read year
    do
      $cdhit/cd-hit -g 1 \
                    -i identification/$geneName\_from_$year\_all.faa \
                    -o dereplication/$geneName\_from_$year.faa \
                    -c 0.97 \
                    -n 5 \
                    -d 0
      $cdhit/clstr2txt.pl dereplication/$geneName\_from_$year.faa.clstr \
        > dereplication/$geneName\_from_$year.tsv
    done
  cd dereplication
  cat $geneName\_from_*.faa > $geneName\_derep.faa
  grep '>' $geneName\_derep.faa | \
    sed 's/>//' \
    > $geneName\_derep_list.txt
  awk -v geneName="$geneName" '{ print $0","geneName }' $geneName\_derep_list.txt >> ../metabolic_gene_key.csv
done


############################################
############################################
# Classify dsrA genes
############################################
############################################

screen -S HCC_dsrA_tree
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=''
scripts=~/HellsCanyon/code/generalUse
# Set up directory
cd ~/HellsCanyon/dataEdited/metabolic_analyses
mkdir sulfur
mkdir sulfur/dsrA

# Copy over my sequences and align them
sed 's/*//' dereplication/dsrA_derep.faa > sulfur/dsrA/dsrA.faa
grep '>' sulfur/dsrA/dsrA.faa | \
  sed 's/>//' \
  > sulfur/dsrA/dsrA_list.txt
metabolic_HMMs=~/HellsCanyon/references/metabolic_HMMs
cd sulfur/dsrA
hmmalign -o dsrA.sto \
            $metabolic_HMMs/TIGR02064.HMM \
            dsrA.faa
$scripts/convert_stockhold_to_fasta.py dsrA.sto

# Copy in reference sequences
cp ~/references/metabolicProteins/sulfur/dsrA/dsrA_karthik_clean.afa \
    dsrA_karthik_clean.afa

# Align seqs
muscle -profile \
        -in1 dsrA_karthik_clean.afa \
        -in2 dsrA.afa \
        -out dsrA_phylogeny.afa
# Trim alignment

# Generate ML tree
#screen -S EG_dsrA_tree
cd ~/HellsCanyon/dataEdited/metabolic_analyses/sulfur/dsrA/
python $scripts/cleanFASTA.py dsrA_phylogeny_masked.afa
mv -f dsrA_phylogeny_masked.afa_temp.fasta dsrA_phylogeny_masked.afa
FastTree dsrA_phylogeny_masked.afa \
    > dsrA_phylogeny_masked.tree



raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s dsrA_phylogeny_masked.afa \
        -n dsrA
