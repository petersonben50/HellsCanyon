#!/usr/bin/env python

#################################
# cleanFASTA.py
# Written by Benjamin D. Peterson

# This script takes the fasta file
# that it is fed and spits out a
# clean version.
# Basically, it removes all line
# breaks internal to the fasta
# sequence, and shortens the header
# fo the first component.
#################################

import os
import sys
import re
from Bio import SeqIO

fasta = sys.argv[1]
output = fasta.rsplit('.fasta')[0] + '_temp.fasta'

# Add the amended part
with open(output, 'w') as outputFile:
    for seq_record in SeqIO.parse(fasta, "fasta"):
        outputFile.write('>' + seq_record.id + '\n')
        outputFile.write(str(seq_record.seq) + '\n')
