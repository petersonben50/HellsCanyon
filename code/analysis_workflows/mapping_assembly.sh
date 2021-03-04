#!/bin/sh

######################
# code/2017_analysis_assembly/assembly_mapping_2017.sh
# Benjamin D. Peterson
######################

############################################
############################################
# Calculate total coverage of each metagenome
############################################
############################################

screen -S HCC_MG_coverage_counting
mkdir ~/HellsCanyon/dataEdited/mapping/
cd ~/HellsCanyon/dataEdited/mapping/
code=~/HellsCanyon/code/generalUse/readfq-master
read_storage=~/HellsCanyon/dataEdited/metagenomes
IFS=$'\n'

echo -e "metagenomeID\tR1\tR2\tsingle\tmerged" > metagenome_coverage.tsv
awk -F ',' '{ print $1 }' ~/HellsCanyon/metadata/lists/metagenome_key.csv | while read metagenome
do
  echo "Counting coverage in" $metagenome
  R1_count=$($code/kseq_fastq_base $read_storage/$metagenome\_R1.fastq.gz | \
                awk -F " " '{ print $5 }')
  R2_count=$($code/kseq_fastq_base $read_storage/$metagenome\_R2.fastq.gz | \
                awk -F " " '{ print $5 }')
  single_count=$($code/kseq_fastq_base $read_storage/$metagenome\_single.fastq.gz | \
                awk -F " " '{ print $5 }')
  merged_count=$($code/kseq_fastq_base $read_storage/$metagenome\_merged.fastq.gz | \
                awk -F " " '{ print $5 }')
  echo -e $metagenome"\t"$R1_count"\t"$R2_count"\t"$single_count"\t"$merged_count >> metagenome_coverage.tsv
done


############################################
############################################
# Mapping reads to scaffolds
############################################
############################################

mkdir ~/HellsCanyon/dataEdited/mapping

######################
# Generate indices
######################
screen -S HCC_assembly_indexing
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PERL5LIB=""
mkdir ~/HellsCanyon/dataEdited/mapping/indices
cd ~/HellsCanyon/dataEdited/mapping/indices
scaffolds=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/assemblies/scaffolds

awk -F ',' '{ print $1 }' ~/HellsCanyon/metadata/lists/assembly_key.csv | while read assembly
do
  if [ -e $scaffolds/$assembly\_assembly.fna ]; then
    echo $assembly "has been cleaned, let's map to it"
    if [ ! -e bowtie_index_$assembly.1.bt2 ]; then
      echo "Building index for" $assembly
      bowtie2-build $scaffolds/$assembly\_assembly.fna \
                    bowtie_index_$assembly
    else
      echo "Already built index for" $assembly
    fi
  else
    echo $assembly "has not been cleaned yet"
  fi
done


######################
# Map reads and process output
######################

# Identify mapping pairs that have not been completed
cd ~/HellsCanyon/dataEdited/mapping/
rm -f ~/HellsCanyon/metadata/lists/mapping_key_pairsRemaining.tsv
touch ~/HellsCanyon/metadata/lists/mapping_key_pairsRemaining.tsv
cat ~/HellsCanyon/metadata/lists/mapping_key.tsv | while read -r line
do
  metagenome=`echo $line | awk -F '\t' '{ print $1 }'`
  assembly=`echo $line | awk -F '\t' '{ print $2 }'`
  if [ -e $metagenome\_to_$assembly.bam ]; then
    echo $metagenome "mapping to" $assembly "already done"
  else
    echo $metagenome "mapping to" $assembly "still needs to be done"
    echo $metagenome", "$assembly >> ~/HellsCanyon/metadata/lists/mapping_key_pairsRemaining.tsv
  fi
done


# Run mapping
chmod +x /home/GLBRCORG/bpeterson26/HellsCanyon/code/executables/read_mapping_assemblies.sh
cd ~/HellsCanyon/dataEdited/mapping/
condor_submit chmod +x /home/GLBRCORG/bpeterson26/HellsCanyon/code/submission/read_mapping_assemblies.sub
