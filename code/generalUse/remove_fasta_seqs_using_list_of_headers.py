#!/usr/bin/env python3.6

################################
# remove_fasta_seqs_using_list_of_headers.py
# Written by Benjamin D. Peterson

# Reads in a fasta sequence and a
# list of headers and returns a
# fasta file without the fasta
# sequences with that header

# usage 'remove_fasta_seqs_using_list_of_headers.py <fasta.file> <list_of_headers> <desired_output_name'
################################

################################
# Load libraries
################################
import os
import sys
from Bio import SeqIO

################################
# Read command line input
################################
fastaFileName = sys.argv[1]
headerFileName = sys.argv[2]
outputFileName = sys.argv[3]

################################
# Pull out needed fasta sequences
################################
with open(outputFileName, "w") as outputFile, open(headerFileName, "r") as remover:
    remove = remover.read().split("\n")
    for seq in SeqIO.parse(fastaFileName, 'fasta'):
        if seq.id not in remove:
            SeqIO.write(seq, outputFile, "fasta")
