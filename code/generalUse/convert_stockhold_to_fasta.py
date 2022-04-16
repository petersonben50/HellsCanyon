#!/usr/bin/env python3.4

################################
# convert_stockholm_fasta.py
# Written by Benjamin D. Peterson

# Reads in an alignment in the Stockholm
# format (with a .sto extension) and
# outputs it in fasta format (.afa).

# usage 'python convert_stockholm_fasta.py <alignment_file.sto>'
################################

import os
import sys
import Bio
from Bio import AlignIO

# Set needed variables
sto = sys.argv[1]
outname = sto.replace('.sto', '.afa')

# Read in fasta alignment
alignment = AlignIO.read(sto, "stockholm")

# Write out Stockholm alignment
AlignIO.write(alignment, outname, "fasta")
