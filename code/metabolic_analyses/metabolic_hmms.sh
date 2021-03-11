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

    awk -F ',' '{ print $1 }' ~/HellsCanyon/metadata/lists/assembly_key.csv | while read assembly
    do
      if [ ! -e identification/$geneName/$assembly.out ]; then
        echo "Searching for" $geneName "in" $assembly
        hmmsearch --tblout identification/$geneName/$assembly\_$geneName.out \
                  --cpu 4 \
                  --cut_tc \
                  $metabolic_HMMs/$HMM \
                  $ORFs/$assembly.faa \
                  > identification/$geneName/$assembly\_$geneName.txt
        lineCount=`wc -l < identification/$geneName/$assembly\_$geneName.out`
        if [ $lineCount -eq 13 ]; then
          echo "No" $gene "hits in" $geneName
        else
          echo "Pulling" $geneName "sequences out of" $assembly
          python $scripts/extract_protein_hitting_HMM.py \
                  identification/$geneName/$assembly\_$geneName.out \
                  $ORFs/$assembly.faa \
                  identification/$geneName/$assembly\_$geneName.faa
        fi
      else
        echo "Search for" $geneName "is already done in" $assembly
      fi
    done

    # Aggregate all sequences and align to HMM

    cd identification

    shopt -s nullglob
    for i in $geneName/*$geneName.faa; do FOUND=$i;break;done

    if [ ! -z $FOUND ]; then
      echo "Concatenating and aligning" $geneName
      cat $geneName/*$geneName.faa \
          > $geneName\_all.faa
      hmmalign -o $geneName.sto \
                  $metabolic_HMMs/$HMM \
                  $geneName\_all.faa
      $scripts/convert_stockhold_to_fasta.py $geneName.sto
      grep '>' $geneName\_all.faa | \
        sed 's/>//' \
        > $geneName\_all_list.txt
    else
      echo "No" $geneName "sequences found at all :("
    fi
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
grep '>' -h identification/*.afa | \
  sed 's/>//' | \
  cut -d"_" -f1,2 | \
  sort | \
  uniq \
  > depth/scaffold_all_list.txt

# Pull out all depths

chmod +x /home/GLBRCORG/bpeterson26/HellsCanyon/code/executables/aggregate_depth_proteins.sh
condor_submit /home/GLBRCORG/bpeterson26/HellsCanyon/code/submission/aggregate_depth_proteins.sub
