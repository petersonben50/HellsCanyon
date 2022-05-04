#!/bin/sh

##################################################
##################################################
# code/binning/PVC_details/collect_PVC_references.sh
# Benjamin D. Peterson
##################################################
##################################################
cd ~/HellsCanyon/dataEdited/binAnalysis
mkdir PVC_details
mkdir PVC_details/genomes
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh



#########################
# Dereplicate bins from this project
#########################
cd ~/HellsCanyon/dataEdited/binAnalysis
mkdir PVC_details/genome_dereplication
mkdir PVC_details/genome_dereplication/genomes
genome_dereplication=~/HellsCanyon/dataEdited/binAnalysis/PVC_details/genome_dereplication
originalBins=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binning/bins_hgcA_keepers/DNA
HomeBio=/home/GLBRCORG/bpeterson26/HellsCanyon/code/HomeBio
# Add manually generated hgcA+ genomes in PVC superphylum
cp $originalBins/anvio_hgcA_0220.fna $genome_dereplication/genomes
cp $originalBins/anvio_hgcA_0261.fna $genome_dereplication/genomes
cp $originalBins/anvio_hgcA_0040.fna $genome_dereplication/genomes
cp $originalBins/anvio_hgcA_0110.fna $genome_dereplication/genomes
# Add autogenerated bins in PVC superphylum
autoBins=/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binning/autoBinning/hqBinSet
grep "Lentisphaeria" $autoBins/taxonomy/taxonomy_summary.txt | awk -F '\t' '{ print $1 }' | while read binID
do
  echo "Copying over" $binID
  cp $autoBins/DNA/$binID.fna $genome_dereplication/genomes/$binID.fna
done
grep "Kiritimatiellae" $autoBins/taxonomy/taxonomy_summary.txt | awk -F '\t' '{ print $1 }' | while read binID
do
  echo "Copying over" $binID
  cp $autoBins/DNA/$binID.fna $genome_dereplication/genomes/$binID.fna
done

cd $genome_dereplication
ls genomes | sed 's/.fna//' > genome_list.txt
mkdir ORFs
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

cat genome_list.txt | while read accessionID
do
  bash $HomeBio/PM4_binGeneration/IMMA_ORF_STAN_BINS.sh -b $accessionID \
                                                        -i genomes/$accessionID.fna \
                                                        -o ORFs \
                                                        -c $HomeBio/fasta_manipulation/cleanFASTA.py
done


# Run ANI comparisons on genomes
# The following workflow is the code to run
# Sarah Stevens's ANI calculator on a
# folder full of bins.
# Details: https://github.com/sstevens2/ani_compare_dag
# Done on submit node
genome_dereplication=~/HellsCanyon/dataEdited/binAnalysis/PVC_details/genome_dereplication
mkdir $genome_dereplication/ANI_comparison
cd $genome_dereplication/ANI_comparison
wget https://ani.jgi-psf.org/download_files/ANIcalculator_v1.tgz
tar -xzvf ANIcalculator_v1.tgz
git clone https://github.com/sstevens2/ani_compare_dag.git
mv ani_compare_dag HCC_PVC_bins
cd $genome_dereplication/ANI_comparison/HCC_PVC_bins/
mkdir bins_to_dereplicate

cp $genome_dereplication/ORFs/*fna bins_to_dereplicate/
echo 'bins_to_dereplicate' > groupslist.txt
# Change path of executable and
# transfer_input_files lines
sed -i "s/executable = group.sh/executable = \/home\/glbrc.org\/bpeterson26\/HellsCanyon\/dataEdited\/binAnalysis\/PVC_details\/genome_dereplication\/ANI_comparison\/HCC_PVC_bins\/group.sh/" group.sub
sed -i "s/transfer_input_files = \/home\/sstevens2\/ANIcalculator_v1\/ANIcalculator,\/home\/sstevens2\/ANIcalculator_v1\/nsimscan,\$(spllist),\$(totransfer)/transfer_input_files = \/home\/glbrc.org\/bpeterson26\/HellsCanyon\/dataEdited\/binAnalysis\/PVC_details\/genome_dereplication\/ANI_comparison\/ANIcalculator_v1\/ANIcalculator,\/home\/glbrc.org\/bpeterson26\/HellsCanyon\/dataEdited\/binAnalysis\/PVC_details\/genome_dereplication\/ANI_comparison\/ANIcalculator_v1\/nsimscan,\$(spllist),\$(totransfer)/" group.sub
#chmod +x /home/glbrc.org/bpeterson26/HellsCanyon/dataEdited/binAnalysis/PVC_details/genome_dereplication/ANI_comparison/HCC_PVC_bins/group.sh
condor_submit_dag runAllANIcompare.dag
# Download output file (bins_to_dereplicate.all.ani.out.cleaned)
# to my computer:
# /Users/benjaminpeterson/Documents/research/HellsCanyon/dataEdited/bins/binAnalysis/PVC_details/genome_dereplication/bins_to_dereplicate.all.ani.out.cleaned

# Upload list of bins to keep: dataEdited/bins/binAnalysis/PVC_details/genome_dereplication/bins_to_keep.txt
# to here: /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/PVC_details/genome_dereplication/bins_to_keep.txt
cd /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/PVC_details
cat genome_dereplication/bins_to_keep.txt | while read binID
do
  cp genome_dereplication/genomes/$binID.fna genomes/
done



#########################
# Collect 5M PVC hgcA+ bins
#########################
cd ~/HellsCanyon/dataEdited/binAnalysis/PVC_details/genomes
grep -E 'KIR_00|LEN_00' ~/5M/dataEdited/binAnalysis/taxonomy/gtdbtk.bac120.summary.tsv | \
  awk -F '\t' '{ print $1 }' | while read binID
do
  echo "copying over" $binID
  cp ~/5M/dataEdited/binAnalysis/bins_processed/$binID.fna .
done


#########################
# Bins from Jones et al, 2019
#########################
cd ~/HellsCanyon/dataEdited/binAnalysis/PVC_details/genomes
grep 'PVC' ~/references/jonesGenomes/bin_names.tsv | \
  awk -F '\t' '{ print $1 }' | \
  while read binID
do
  echo "cleaning and moving" $binID
  python $HomeBio/fasta_manipulation/cleanFASTA.py ~/references/jonesGenomes/$binID.fna
  mv ~/references/jonesGenomes/$binID.fna_tempCleanedFile ./$binID.fna
done


#########################
# References selected for 5M project
#########################
cd ~/HellsCanyon/dataEdited/binAnalysis/PVC_details/genomes
grep -E 'c__Kiritimatiellae|c__Lentisphaeria|staleyi|formosa' ~/5M/dataEdited/binAnalysis/phylogeny/PVC/reference_taxonomy.tsv | \
  awk -F '\t' '{ print $1 }' | \
  sed 's/GB_//' | \
  sed 's/RS_//' | \
  while read binID
do
  ls ~/references/genomes/bins/$binID*.fna
  cp ~/references/genomes/bins/$binID*.fna ./$binID.fna
  python $HomeBio/fasta_manipulation/cleanFASTA.py $binID.fna
  mv -f $binID.fna_tempCleanedFile $binID.fna
done



#########################
# Additional genomes
#########################
cd ~/HellsCanyon/dataEdited/binAnalysis/PVC_details/genomes
binID=GCF_003096415.1
echo "moving" $binID
cp ~/references/genomes/bins/$binID*.fna ./$binID.fna
python $HomeBio/fasta_manipulation/cleanFASTA.py $binID.fna
mv -f $binID.fna_tempCleanedFile $binID.fna

binID=GCF_900890425.1
echo "moving" $binID
cp ~/references/genomes/bins/$binID*.fna ./$binID.fna
python $HomeBio/fasta_manipulation/cleanFASTA.py $binID.fna
mv -f $binID.fna_tempCleanedFile $binID.fna

binID=GCF_900890705.1
echo "moving" $binID
cp ~/references/genomes/bins/$binID*.fna ./$binID.fna
python $HomeBio/fasta_manipulation/cleanFASTA.py $binID.fna
mv -f $binID.fna_tempCleanedFile $binID.fna


#########################
# Predict ORFs
#########################
cd ~/HellsCanyon/dataEdited/binAnalysis/PVC_details
mkdir ORFs
ls genomes | sed 's/.fna//' > genome_list.txt
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
cat genome_list.txt | while read binID
do
  if [ ! -e ORFs/$binID.faa ]; then
    echo "Working on ORF prediction for" $binID
    bash $HomeBio/PM4_binGeneration/IMMA_ORF_STAN_BINS.sh -b $binID \
                                                          -i genomes/$binID.fna \
                                                          -o ORFs \
                                                          -c $HomeBio/fasta_manipulation/cleanFASTA.py
  else
    echo "ORF prediction for" $binID "already done"
  fi
done
