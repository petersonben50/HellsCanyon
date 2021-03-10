#!/bin/sh

########################
# code/mer/test_pfam_hmm.sh
# Benjamin D. Peterson

# This script will be used to
# download and test the pfam
# HMM for merB.
########################

# Check for cut-offs
head -n 30 /Users/benjaminpeterson/Documents/research/HellsCanyon/references/merB/MerB.hmm


# Test out HMM

conda activate bioinformatics
HCC
mkdir dataEdited/mer/exploration/HMM_test
hmmsearch --tblout dataEdited/mer/exploration/HMM_test/merB.out \
          --cpu 4 \
          --cut_tc \
          references/merB/MerB.hmm \
          dataEdited/mer/exploration/uniprot_ec_4.99.1.2.fasta \
          > dataEdited/mer/exploration/HMM_test/merB_output.txt

hmmsearch --tblout dataEdited/mer/exploration/HMM_test/merB_noCutoff.out \
          --cpu 4 \
          -T 1 \
          references/merB/MerB.hmm \
          dataEdited/mer/exploration/uniprot_ec_4.99.1.2.fasta \
          > dataEdited/mer/exploration/HMM_test/merB_noCutoff_output.txt


########################
# Pull out 10 seqs that don't hit HMM
########################
HCC
grep '>' dataEdited/mer/exploration/uniprot_ec_4.99.1.2.fasta | \
  sed 's/>//' | \
  awk '{ print $1 }' \
  > dataEdited/mer/exploration/HMM_test/proteins_all.txt

cd dataEdited/mer/exploration/HMM_test/
grep -v "#" merB_noCutoff.out | \
  awk '{ print $1 }' \
  > proteins_that_hit_hmm.txt

cat proteins_all.txt proteins_that_hit_hmm.txt | \
  sort | uniq -u
