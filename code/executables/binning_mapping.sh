#!/bin/sh

#########################
# binning_mapping.sh
# Benjamin D. Peterson

# This script maps the reads from each
# metagenome to the assemblies from this analysis.
# It includes the index building step.
# I used bowtie2 (v2.2.2) for this.
# The submission file loops over each assembly.
#########################


# Will read in these variables from
# submission file:
# $assembly

# Will read in this path from submission file:
# mapping=/home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_binning/binning_initial/mapping


#########################
# Activate environment
#########################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
cd $mapping


#########################
# Build indicies
#########################
if [ ! -e $mapping/indices/$assembly.1.bt2 ]; then
  echo "Make index for assembly" $assembly
  /opt/bifxapps/bowtie2-2.2.2/bowtie2-build $scaffoldsLocation/$assembly\_filtered_scaffolds.fna \
                                            $mapping/indices/$assembly
else
  echo "Already made index for" $assembly
fi


#########################
# Set up list of metagenomes for mapping
#########################
awk -F '\t' -v assembly="$assembly" '$2 == assembly { print $1 }' $mappingKey > metagenomes_for_$assembly\_mapping.txt


#########################
# Do mapping
#########################
IFS=$'\n'
cat metagenomes_for_$assembly\_mapping.txt | while read metagenome
do
  if [ ! -e $mapping/$metagenome\_to_$assembly.bam ]; then
    echo "Mapping MG" $metagenome "to" $assembly
    /opt/bifxapps/bowtie2-2.2.2/bowtie2 -x $mapping/indices/$assembly \
                                        -1 metagenomes/$metagenome\_R1.fastq.gz \
                                        -2 metagenomes/$metagenome\_R2.fastq.gz \
                                        -U metagenomes/$metagenome\_single.fastq.gz,metagenomes/$metagenome\_merged.fastq.gz \
                                        -p 10 \
                                        -S $mapping/$metagenome\_to_$assembly.sam
    echo "Converting, sorting, and indexing mapping data from" $metagenome "to" $assembly
    samtools view $mapping/$metagenome\_to_$assembly.sam \
                  -o $mapping/raw.$metagenome\_to_$assembly.bam
    samtools sort -m 10G \
                  -@ 10 \
                  $mapping/raw.$metagenome\_to_$assembly.bam \
                  -o $mapping/$metagenome\_to_$assembly.bam
    samtools index $mapping/$metagenome\_to_$assembly.bam

    rm -f $mapping/$metagenome\_to_$assembly.sam \
          $mapping/raw.$metagenome\_to_$assembly.bam

  else
    echo "Already mapped" $metagenome "to" $assembly
  fi
done

rm -f metagenomes_for_$assembly\_mapping.txt
