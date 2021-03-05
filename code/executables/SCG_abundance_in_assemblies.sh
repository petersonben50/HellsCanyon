#!/bin/sh

#########################
# SCG_abundance_in_assemblies.sh
# Benjamin D. Peterson

# This script calculates the read coverage
# from a series of metagenomes over a set of
# single copy genes
#########################


# Will read in these variables from
# submission file:
# geneName,scripts,scgHmms_location,scgHmms_key,assembly_list
# assembly_location,metagenome_list,metagenome_location


#########################
# Activate environment
#########################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''


#########################
# Set up directory
#########################
cd $output_directory
if [ -d working_directory_$geneName ]; then
  rm -rf working_directory_$geneName
fi
mkdir working_directory_$geneName
cd working_directory_$geneName


#########################
# Pull out HMM names
#########################
# For some SCGs, there might be different bacterial
# and archaeal HMMs. So, we'll search the assemblies
# using both and combine the abundances.
awk -F ',' -v geneName="$geneName" ' $1 == geneName { print $2 }' $scgHmms_key \
  > hmm_list.txt


#########################
# Search for SCGs in all assemblies
#########################
cat $assembly_list | while read assembly
do
  cat hmm_list.txt | while read hmmName
  do
    echo "Searching" $assembly "with" $hmmName
    hmmsearch --tblout $assembly\_$geneName\_$hmmName.out \
              --cpu 6 \
              --cut_nc \
              $scgHmms_location/$hmmName \
              $assembly_location/$assembly.faa \
              > $assembly\_$geneName\_$hmmName.txt
    lineCount=`wc -l < $assembly\_$geneName\_$hmmName.out`
    if [ $lineCount -eq 13 ]; then
      echo "No" $gene "hits in" $geneName
    else
      echo "Pulling" $geneName "sequences out of" $assembly
      python $scripts/extract_protein_hitting_HMM.py \
              $assembly\_$geneName\_$hmmName.out \
              $assembly_location/$assembly.faa \
              $assembly\_$geneName\_$hmmName.faa
    fi
  done
done


#########################
# Dereplicate SCGs across assemblies within years
#########################
rm -f *_raw.faa
cat $assembly_key | while read line
do
  assembly=`echo $line | awk -F ',' '{ print $1 }'`
  year=`echo $line | awk -F ',' '{ print $2 }'`
  echo $assembly "is from" $year
  cat $assembly\_$geneName\_*.faa >> $geneName\_$year\_raw.faa
done

# Dereplicate sequences within year
rm -f $geneName.faa
cdhit=~/programs/cdhit-master
awk -F ',' '{ print $2 }' $assembly_key | \
  sort | uniq | \
  while read year
  do
    echo "Dereplicating" $geneName "sequences from" $year
    ls $geneName\_$year\_raw.faa
    $cdhit/cd-hit -g 1 \
                  -i $geneName\_$year\_raw.faa \
                  -o $geneName\_$year.faa \
                  -c 0.97 \
                  -n 5 \
                  -d 0
    cat $geneName\_$year.faa >> $geneName.faa
  done

grep '>' $geneName.faa | \
  sed 's/>//' \
  > $geneName\_list.txt
# Generate list of all the scaffolds
cat $geneName\_list.txt | \
  cut -d"_" -f1,2 \
  > $geneName\_scaffold_list.txt


#########################
# Caculate coverage of scaffold with $geneName
#########################
rm -f $geneName\_scg_coverage.tsv
cat $metagenome_list | while read metagenome
do
  conda activate bioinformatics
  PERL5LIB=""
  rm -f $metagenome\_$geneName\_depth_raw.tsv
  # Calculate depth over each residue in each scaffold
  echo "Calculating coverage of each residue in a" $geneName"+ scaffold in" $metagenome

  cat $geneName\_scaffold_list.txt | while read scaffold
  do
    assembly=$(echo $scaffold | awk -F '_' '{ print $1 }')
    if [ -e $metagenome_location/$metagenome\_to_$assembly.bam ]; then
      samtools depth -a -r $scaffold $metagenome_location/$metagenome\_to_$assembly.bam \
          >> $metagenome\_$geneName\_depth_raw.tsv
    else
      echo $metagenome "was not mapped to" $assembly
    fi
  done
  conda deactivate

  # Average the depth of each residue over the entire
  echo "Aggregating coverage of each" $geneName"+ scaffold in" $metagenome
  conda activate py_viz
  PYTHONPATH=""
  python $scripts/calculate_depth_contigs.py \
            $metagenome\_$geneName\_depth_raw.tsv \
            150 \
            $metagenome\_$geneName\_depth.tsv
  conda deactivate
  rm -f $metagenome\_$geneName\_depth_raw.tsv

  # Sum the coverage of all the identified genes in each metagenome
  echo "Summing" $geneName "coverage in" $metagenome
  coverage_count=`awk -F '\t' '{sum += $2} END {print sum}' $metagenome\_$geneName\_depth.tsv`
  echo -e "$geneName\t$metagenome\t$coverage_count" >> $output_directory/$geneName\_scg_coverage.tsv
done
