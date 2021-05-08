#!/bin/sh

######################
# code/assembly_processing.sh
# Benjamin D. Peterson

# This set of scripts will contain everything
# needed for the initial processing of the
# metagenomes from the 2019 intensive.
######################


############################################
############################################
# Data retrieval and inspection
############################################
############################################

######################
# KMBP004: First round of metagenomes
######################

# Get set up
mkdir ~/HellsCanyon/dataRaw/metagenomes/2019_intensive


# Now retrieve the metagenomes using lftp
# (in the lftp virtual environment)
screen -S HCC_metagenome_retrieval
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate lftp

# Download first attempt at sequencing for KMBP004.
#cd ~/HellsCanyon/dataRaw/metagenomes/2019_intensive
#lftp -c 'set ssl:verify-certificate no set ftp:ssl-protect-data true set ftp:ssl-force true; open -u n191204_Peterson,goofu5Dae7voo7f -e "mirror -c; quit" ftp://gslanalyzer.qb3.berkeley.edu:990'
#mv laneBarcode.html KMBP004_laneBarcode.html
#mv md5sum.txt KMBP004_md5sum.txt

# Move crappy metagenomes to own folder
#cd ~/HellsCanyon/dataRaw/metagenomes/2019_intensive
#mkdir misrun_samples
#mv KMBP004[A-D]* misrun_samples

# Download samples in KMBP004 that were re-run.
#cd ~/HellsCanyon/dataRaw/metagenomes/2019_intensive
#mkdir temp
#cd temp
#lftp -c 'set ssl:verify-certificate no set ftp:ssl-protect-data true set ftp:ssl-force true; open -u n191218_Peterson,paiziX3pheeb1ie -e "mirror -c; quit" ftp://gslanalyzer.qb3.berkeley.edu:990'
#mv laneBarcode.html KMBP004_redo_laneBarcode.html
#mv md5sum.txt KMBP004_redo_md5sum.txt
#mv KMBP004* ~/HellsCanyon/dataRaw/metagenomes/2019_intensive
#cd ~/HellsCanyon/dataRaw/metagenomes/2019_intensive
#rm -f -r temp



# Let's get a list of our metagenomes.
# Save out the 2019 intensive metagenome IDs here:
# ~/HellsCanyon/dataEdited/metagenomes/2019_intensive_MG_list.txt
# Save out all metagenome IDs here:
# ~/HellsCanyon/dataEdited/metagenomes/MG_list_all.txt



######################
# KMBP007, KMBP008, KMBP009: Summer 2020 sequencing submissions
######################

# Log onto GLBRC server
#screen -S HCC_MG_transfer
#cd ~/HellsCanyon/dataRaw/metagenomes

#source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
#conda activate bioinformatics
#wget -r -np -nH --ask-password ftps://mcmahonk@titan.bch.msu.edu/20200731_DNASeq_PE150
#mv 20200731_DNASeq_PE150/* ./
#rm -fr 20200731_DNASeq_PE150
#exit

############################################
############################################
# Trimming metagenomes
############################################
############################################

screen -S HCC_metagenome_trimming
mkdir ~/HellsCanyon/dataEdited/metagenomes
mkdir ~/HellsCanyon/dataEdited/metagenomes/reports

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PERL5LIB=""

######################
# Trim metagenomes
######################

cd ~/HellsCanyon/dataRaw/metagenomes

read_storage=~/HellsCanyon/dataEdited/metagenomes
ancillary_info=~/HellsCanyon/dataEdited/metagenomes/reports

cat ~/HellsCanyon/metadata/lists/MG_list_all.txt | while read metagenome
do
  if [ ! -e $read_storage/$metagenome\_R1.fastq.gz ]; then
    echo "Processing" $metagenome
    fastp --in1 $metagenome*R1*fastq.gz \
          --in2 $metagenome*R2*fastq.gz \
          --out1 $read_storage/$metagenome\_R1.fastq.gz \
          --out2 $read_storage/$metagenome\_R2.fastq.gz \
          --unpaired1 $read_storage/$metagenome\_single.fastq.gz \
          --unpaired2 $read_storage/$metagenome\_single.fastq.gz \
          --merge \
          --merged_out $read_storage/$metagenome\_merged.fastq.gz \
          --failed_out $ancillary_info/$metagenome\_failed.fastq.gz \
          --html $ancillary_info/$metagenome\_report.html \
          --detect_adapter_for_pe \
          --cut_tail \
          --cut_tail_window_size 10 \
          --cut_tail_mean_quality 20 \
          --length_required 100
    fastqc -o $ancillary_info \
            $read_storage/$metagenome\_R1.fastq.gz
    fastqc -o $ancillary_info \
            $read_storage/$metagenome\_R2.fastq.gz
    fastqc -o $ancillary_info \
            $read_storage/$metagenome\_single.fastq.gz

  else
    echo "Already processed" $metagenome
  fi
done



######################
# Count reads in metagenome pre-trimming
######################

screen -S HCC_metagenome_read_counting
cd ~/HellsCanyon/dataRaw/metagenomes
read_storage=~/HellsCanyon/dataRaw/metagenomes
ancillary_info=~/HellsCanyon/dataEdited/metagenomes/reports

echo -e "metagenomeID\tforwardReads\treverseReads" > $ancillary_info/metagenome_read_count_pretrim.tsv
cat ~/HellsCanyon/metadata/lists/MG_list_all.txt | while read metagenome
do
  echo "Working on" $metagenome
  forwardCount=$(zgrep -c "@" $read_storage/$metagenome*_R1_001.fastq.gz)
  reverseCount=$(zgrep -c "@" $read_storage/$metagenome*_R2_001.fastq.gz)
  echo -e $metagenome"\t"$forwardCount"\t"$reverseCount >> $ancillary_info/metagenome_read_count_pretrim.tsv
done






######################
# Count reads in metagenome
######################

screen -S HCC_metagenome_read_counting
cd ~/HellsCanyon/dataRaw/metagenomes
read_storage=~/HellsCanyon/dataEdited/metagenomes
ancillary_info=~/HellsCanyon/dataEdited/metagenomes/reports

echo -e "metagenomeID\tforwardReads\treverseReads\tsingleReads\tmergedReads" > $ancillary_info/metagenome_read_count.tsv
cat ~/HellsCanyon/metadata/lists/MG_list_all.txt | while read metagenome
do
  echo "Working on" $metagenome
  forwardCount=$(zgrep -c "@" $read_storage/$metagenome\_R1.fastq.gz)
  reverseCount=$(zgrep -c "@" $read_storage/$metagenome\_R2.fastq.gz)
  singleCount=$(zgrep -c "@" $read_storage/$metagenome\_single.fastq.gz)
  mergeCount=$(zgrep -c "@" $read_storage/$metagenome\_merged.fastq.gz)
  echo -e $metagenome"\t"$forwardCount"\t"$reverseCount"\t"$singleCount"\t"$mergeCount >> $ancillary_info/metagenome_read_count.tsv
done








######################
# Coverage post-trimming
######################

screen -S HCC_MG_coverage_counting_post
cd ~/HellsCanyon/dataEdited/metagenomes/reports
code=~/HellsCanyon/code/generalUse/readfq-master
read_storage=~/HellsCanyon/dataEdited/metagenomes
ancillary_info=~/HellsCanyon/dataEdited/metagenomes/reports

echo -e "metagenomeID\tR1\tR2\tsingle\tmerged" > $ancillary_info/metagenome_coverage.tsv
cat ~/HellsCanyon/metadata/lists/MG_list_all.txt | while read metagenome
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
  echo -e $metagenome"\t"$R1_count"\t"$R2_count"\t"$single_count"\t"$merged_count >> $ancillary_info/metagenome_coverage.tsv
done









######################
# Cluster metagenomes using Mash
######################

screen -S HCC_mash
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics


# Generate sketches for each metagenome

mkdir ~/HellsCanyon/dataEdited/mash_data
mkdir ~/HellsCanyon/dataEdited/mash_data/temp_MG_files
mkdir ~/HellsCanyon/dataEdited/mash_data/sketch_files
cd ~/HellsCanyon/dataEdited/mash_data
read_storage=~/HellsCanyon/dataEdited/metagenomes


cat ~/HellsCanyon/metadata/lists/MG_list_all.txt | while read metagenome
do
  if [ ! -e sketch_files/$metagenome.msh ]; then
    cat $read_storage/$metagenome*.fastq.gz > temp_MG_files/$metagenome.fastq.gz
    mash sketch -S 50 \
                -r \
                -m 2 \
                -k 21 \
                -s 100000 \
                -o sketch_files/$metagenome \
                temp_MG_files/$metagenome.fastq.gz
    rm -f temp_MG_files/$metagenome.fastq.gz
  else
    echo "Already sketched" $metagenome
  fi
done



# Merge sketches

cd ~/HellsCanyon/dataEdited/mash_data
mash paste sketch_files/HCC_MG_all sketch_files/KMBP00*.msh

mash paste sketch_files/HCC_MG_2017 \
            sketch_files/KMBP002*.msh \
            sketch_files/KMBP003*.msh

mash paste sketch_files/HCC_MG_2018 \
            sketch_files/KMBP007*.msh \
            sketch_files/KMBP008*.msh

mash paste sketch_files/HCC_MG_2019 \
            sketch_files/KMBP004*.msh \
            sketch_files/KMBP009*.msh

# Run distance calcuation

mash dist -S 50 \
          sketch_files/HCC_MG_all.msh \
          sketch_files/HCC_MG_all.msh \
          > HCC_MG_all.dist

mash dist -S 50 \
          sketch_files/HCC_MG_2017.msh \
          sketch_files/HCC_MG_2017.msh \
          > HCC_MG_2017.dist

mash dist -S 50 \
          sketch_files/HCC_MG_2018.msh \
          sketch_files/HCC_MG_2018.msh \
          > HCC_MG_2018.dist

mash dist -S 50 \
          sketch_files/HCC_MG_2019.msh \
          sketch_files/HCC_MG_2019.msh \
          > HCC_MG_2019.dist




############################################
############################################
# Metagenome assembly
############################################
############################################



######################
# Assemblies by metaSPades
######################

# Make sure the list of wanted assemblies by metaSPADes is up to date
# in ~/HellsCanyon/metadata/lists/assembly_group_metaspades.csv

# Get a list of group names
cd  ~/HellsCanyon/metadata/lists
tail -n +2 assembly_group_metaspades.csv | \
    awk -F ',' '{ print $1 }' | \
    sort | uniq \
    > assembly_list_metaspades.txt



screen -S HCC_metagenome_coassembly

cd ~/HellsCanyon/dataEdited/assemblies

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""

code=~/HellsCanyon/code/assembly_processing/
assembly_grouping=~/HellsCanyon/metadata/lists/assembly_group_metaspades.csv
read_storage=~/HellsCanyon/dataEdited/metagenomes/
output=~/HellsCanyon/dataEdited/assemblies/assembly_files/


cat ~/HellsCanyon/metadata/lists/assembly_list_metaspades.txt | while read assembly
do

  if [ ! -d $output/$assembly ]; then
    mkdir $output/$assembly
  fi

  if [ ! -e $output/$assembly/scaffolds.fasta ]; then
    echo "Assembling" $assembly
    python $code/assembly_by_group_metaspades.py $assembly \
                                                  $assembly_grouping \
                                                  $read_storage \
                                                  $output/$assembly
    # To continue a paused run:
    # assembly=july2019_highredox_metaspades
    # metaspades.py --continue -o assembly

  else
    echo $assembly "already assembled"
  fi
done


#exit

######################
# Clean up assemblies for metaSPADes
######################

screen -S HCC_clean_metagenome_assembly

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=/home/GLBRCORG/bpeterson26/miniconda3/envs/anvio5/lib/python3.6/site-packages/

cd ~/HellsCanyon/dataEdited/assemblies
mkdir scaffolds
mkdir scaffolds/renaming_reports


cat ~/HellsCanyon/metadata/lists/assembly_list_metaspades.txt | while read assembly
do
  if [ -e assembly_files/$assembly/scaffolds.fasta ]; then
    echo $assembly "is done assembling. Let's clean it."
    if [ -e scaffolds/$assembly\_assembly.fna ]; then
      echo "Dude, you already cleaned the assembly for" $assembly". Relax."
    else
      echo "Cleaning the assembly for" $assembly
      anvi-script-reformat-fasta assembly_files/$assembly/scaffolds.fasta \
                                  -o scaffolds/$assembly\_assembly.fna \
                                  -l 1000 \
                                  --simplify-names \
                                  --prefix $assembly \
                                  --report-file scaffolds/renaming_reports/$assembly\_report_file.txt
    fi
  else
    echo "Yo, you gotta go back and assemble" $assembly
  fi
done

conda deactivate





######################
# Combine lists of assemblies
######################
cd ~/HellsCanyon/metadata/lists
cat assembly_list_megahit.txt \
    assembly_list_metaspades.txt \
    > assembly_list.txt


#exit







######################
# Get assembly stats
######################

screen -S HCC_assembly_stats

cd ~/HellsCanyon/dataEdited/assemblies
mkdir reports
scripts=~/HellsCanyon/code/assembly_processing

cat ~/HellsCanyon/metadata/lists/assembly_list.txt | while read assembly
do
  if [ -e scaffolds/$assembly\_assembly.fna ]; then
    echo $assembly "has been cleaned. Let's get some stats on it."
    if [ ! -e reports/$assembly\_report.txt ]; then
      perl $scripts/abyss-fac.pl scaffolds/$assembly\_assembly.fna \
          > reports/$assembly\_report.txt
    else
      echo "Already checked" $assembly
    fi
  else
    echo $assembly "has not been cleaned. Do that first."
  fi
done

# Aggregate stats
cd ~/HellsCanyon/dataEdited/assemblies/reports
echo -e "n\tn:200\tL50\tmin\tN80\tN50\tN20\tmax\tsum\tassemblyID" > all_assemblies_stats.txt
for file in *report.txt
do
  tail -n +2 $file >> all_assemblies_stats.txt
done
# rm KMBP*_report.txt

# Clean up the report file
cat ~/HellsCanyon/metadata/lists/assembly_list.txt | while read assembly
do
  sed "s/scaffolds\/$assembly\_assembly.fna/$assembly/" all_assemblies_stats.txt \
    > all_assemblies_stats.txt_edited
  mv -f all_assemblies_stats.txt_edited all_assemblies_stats.txt
done

# exit








############################################
############################################
# Predict ORFs
############################################
############################################

mkdir /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/assemblies/ORFs
chmod +x /home/GLBRCORG/bpeterson26/HellsCanyon/code/executables/ORF_prediction_assemblies.sh
condor_submit /home/GLBRCORG/bpeterson26/HellsCanyon/code/submission/ORF_prediction_assemblies.sub











######################
# Count ORFs
######################

cd ~/HellsCanyon/dataEdited/assemblies

echo -e 'assemblyID\tORF_count' > reports/ORF_counts.tsv
cat ~/HellsCanyon/metadata/lists/assembly_list.txt | while read assembly
do
  if [ -e ORFs/$assembly.faa ]; then
    echo "Count ORFs in" $assembly
    ORF_count=$(grep '>' ORFs/$assembly.faa | wc -l)
    echo -e "$assembly\t$ORF_count" >> reports/ORF_counts.tsv
  fi
done
