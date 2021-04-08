#################################
# batch_HMMs.py
# Benjamin D. Peterson

# This python script will run a set
# of HMMs, as defined by a csv file,
# against a set of ORFs concatentated
# from a set of bins.
#################################

# Set up environment
import os
import sys
import argparse
import glob
from Bio import AlignIO
from Bio import SearchIO
from Bio import SeqIO
import pandas as pd
import subprocess

# Set up an argument parser
parser = argparse.ArgumentParser()
parser.add_argument('--orf_file')
parser.add_argument('--bin_folder')
parser.add_argument('--orf_folder')
parser.add_argument('--g2b')
parser.add_argument('--hmm_folder')
parser.add_argument('--hmm_csv')
parser.add_argument('--output')
parser.add_argument('--testing')

# Parse names from argument
inputs = parser.parse_args()
ORF_FILE = inputs.orf_file
BIN_FOLDER = inputs.bin_folder
ORF_FOLDER = inputs.orf_folder
G2B = inputs.g2b
HMM_FOLDER = inputs.hmm_folder
HMM_CSV = inputs.hmm_csv
OUTPUT_LOCATION = inputs.output
TESTING = inputs.testing


###########################
# Make sure we're set up properly
###########################

# Make sure we have either an ORF folder or ORF file + G2B file
if ORF_FILE is not None:
    if ORF_FOLDER is not None:
        print("You've supplied both an concatenated ORF file and a folder with the ORFs for the bins. Only supply one.")
        sys.exit()
    elif BIN_FOLDER is not None:
        print("You've supplied both an concatenated ORF file and a folder with the bins. Only supply one.")
        sys.exit()
    elif G2B is None:
        print("You supplied a concatenated ORF file, you'll also want to supply a gene-to-bin file.")
        sys.exit()
    else:
        print("Input is concatenated ORFs file")
        input_type = "concat_orfs"
elif ORF_FOLDER is not None:
    if BIN_FOLDER is not None:
        print("You've supplied both a folder with the bins and a folder with ORFs. Only supply one.")
        sys.exit()
    else:
        print("Input is a folder of ORF files")
        input_type = "folder_of_orfs"
elif BIN_FOLDER is not None:
    print("Input is a folder of bin files")
    input_type = "folder_of_bins"
else:
    print("You haven't supplied an input for the bins")
    sys.exit()

# Check if output directory already exists
if os.path.isdir(OUTPUT_LOCATION) == True:
    print("Hey dummy, " + OUTPUT_LOCATION + " is already a directory. You're gonna write over stuff.")
    sys.exit()
os.mkdir(OUTPUT_LOCATION)


###########################
# Get inputs set up
###########################
# Set up inputs
FAA_INPUTS = OUTPUT_LOCATION + '/faa_inputs/'
concat_orf_to_use = FAA_INPUTS + 'ORFs.faa'
os.mkdir(FAA_INPUTS)


###########################
# Input set-up: using folder of bins
###########################
ORF_FILES = OUTPUT_LOCATION + '/ORF_files/'
os.mkdir(ORF_FILES)
if input_type == "folder_of_bins":
    print("Running ORF prediction for bins")
    genome_files = BIN_FOLDER + "/**"
    genomes = glob.glob(genome_files)
    for genome in genomes:
        if genome.endswith('.fna'):
            genome_name = genome.rsplit('.', 1)[0].rsplit('/', 1)[1]
            prodigal_cmd = 'prodigal -i ' + genome
            prodigal_cmd = prodigal_cmd + ' -o ' + ORF_FILES + genome_name + '.gff -f gff '
            prodigal_cmd = prodigal_cmd + '-a ' + ORF_FILES + genome_name + '.faa '
            prodigal_cmd = prodigal_cmd + ' -p single -q'
            os.system(prodigal_cmd)
    input_type = "folder_of_orfs"
    ORF_FOLDER = ORF_FILES


###########################
# Input set-up: using folder of ORFs
###########################
# If we supplied folder of ORFs, concatenate them to a folder, and generate G2B file
if input_type == "folder_of_orfs":
    print("Concatenating ORFs from all bins")
    concat_cmd = "cat "
    genome_files = ORF_FOLDER + "/**"
    genomes = glob.glob(genome_files)
    for genome in genomes:
        if genome.endswith('.faa'):
            concat_cmd = concat_cmd + " " + genome
    concat_cmd = concat_cmd + " > " + concat_orf_to_use
    os.system(concat_cmd)

# Generate gene to bin file
if input_type == "folder_of_orfs":
    print("Generating G2B file from all bins")
    g2bkey = dict()
    genome_files = ORF_FOLDER + "/**"
    genomes = glob.glob(genome_files)
    for genome in genomes:
        if genome.endswith('.faa'):
            genome_name = genome.rsplit("/", 1)[1].rsplit(".faa", 1)[0]
            genome_orfs = SeqIO.parse(genome, "fasta")
            for ORF in genome_orfs:
                g2bkey[ORF.id] = genome_name


###########################
# Input set-up: using concatenated ORF file
###########################
# If we supplied concatenated ORF, move it over:
if input_type == "concat_orfs":
    print("Copying " + ORF_FILE + " to " + concat_orf_to_use)
    cp_command = "cp " + ORF_FILE + " " + concat_orf_to_use
    os.system(cp_command)

# Set up G2B file
if input_type == "concat_orfs":
    print("Setting up G2B dictionary")
    g2b = pd.read_csv(G2B, delimiter="\t", names=['gene', 'bin'])
    g2bkey = dict()
    for index, row in g2b.iterrows():
        g2bkey[row.gene] = row.bin


###########################
# Read in needed files
###########################
# Read in dataframe of HMM info
hmm_csv = pd.read_csv(HMM_CSV, delimiter=",")
# Set up dictionary for fasta sequences
faadict = dict()
for seq_record in SeqIO.parse(concat_orf_to_use, "fasta"):
	faadict[seq_record.id] = seq_record.seq


###########################
# Set up output directories
###########################
HMM_OUTPUT = OUTPUT_LOCATION + "/hmm_output/"
os.mkdir(HMM_OUTPUT)
ALIGNMENT_OUTPUT = OUTPUT_LOCATION + "/alignments/"
os.mkdir(ALIGNMENT_OUTPUT)
HMM_HITS = OUTPUT_LOCATION + "/alignments/hmm_hits/"
os.mkdir(HMM_HITS)
BIN_COUNTS = OUTPUT_LOCATION + "/bin_counts/"
os.mkdir(BIN_COUNTS)


###########################
# Lines to break script while testing, if needed
###########################
if TESTING:
    print("Oh " + TESTING + " I'm testing")
    sys.exit()


###########################
# Run all HMMs on ORF file
###########################
for index, row in hmm_csv.iterrows():
    prot_name = row['protein']
    hmm_id = row['hmmID']
    HMM = HMM_FOLDER + "/" + hmm_id
    print("Searching for " + prot_name)
    hmm_cmd = 'hmmsearch --tblout ' + HMM_OUTPUT + prot_name + '.out --cpu 4 --cut_tc ' + HMM + " " + concat_orf_to_use + " > " + HMM_OUTPUT + prot_name + '.txt'
    os.system(hmm_cmd)


###########################
# Extract amino acid sequences of all hits to each HMMs
###########################
for index, row in hmm_csv.iterrows():
    prot_name = row['protein']
    HMMER_OUTPUT_FILE = HMM_OUTPUT + prot_name + '.out'
    HMM_OUTPUT_LENGTH = subprocess.check_output('wc -l < ' + HMMER_OUTPUT_FILE, shell=True)
    if int(HMM_OUTPUT_LENGTH) > 13:
        hmmer_output = SearchIO.read(HMMER_OUTPUT_FILE, 'hmmer3-tab')
        fasta_output_for_hits = HMM_HITS + prot_name + ".faa"
        print("Extracting AA sequences for " + prot_name)
        with open(fasta_output_for_hits, 'w') as resultFile:
            for sampleID in hmmer_output:
                resultFile.write('>' + str(sampleID.id) + ' ' + str(sampleID.bitscore) + '\n' + str(faadict[sampleID.id]) + '\n')
    else:
        print('No hits for ' + prot_name)


###########################
# Align amino acid sequences to HMM
###########################
for index, row in hmm_csv.iterrows():
    prot_name = row['protein']
    hmm_id = row['hmmID']
    sto_output = ALIGNMENT_OUTPUT + prot_name + '.sto'
    afa_output = ALIGNMENT_OUTPUT + prot_name + '.afa'
    HMM = HMM_FOLDER + "/" + hmm_id
    fasta_output_for_hits = HMM_HITS + prot_name + ".faa"
    if os.path.isfile(fasta_output_for_hits):
        print("Aligning sequences of " + prot_name + " to HMM")
        hmmalign_cmd = 'hmmalign -o ' + sto_output + ' ' + HMM + " " + fasta_output_for_hits
        os.system(hmmalign_cmd)
        # Read in fasta alignment
        sto_alignment = AlignIO.read(sto_output, "stockholm")
        # Write out afa alignment
        AlignIO.write(sto_alignment, afa_output, "fasta")
        os.remove(sto_output)


###########################
# Generate dataframe with hits against bins
###########################
for index, row in hmm_csv.iterrows():
    bin_hits = dict()
    prot_name = row['protein']
    HMMER_OUTPUT_FILE = HMM_OUTPUT + prot_name + '.out'
    HMM_OUTPUT_LENGTH = subprocess.check_output('wc -l < ' + HMMER_OUTPUT_FILE, shell=True)
    if int(HMM_OUTPUT_LENGTH) > 13:
        hmmer_output = SearchIO.read(HMMER_OUTPUT_FILE, 'hmmer3-tab')
        print("Counting bin hits for " + prot_name)
        for sampleID in hmmer_output:
            bin_hits[sampleID.id] = [sampleID.id, g2bkey[sampleID.id], prot_name]
            bin_hits_df = pd.DataFrame.from_dict(bin_hits, orient = 'index', columns = ['geneID', 'binID', "proteinName"])
        if 'all_bin_hits' in locals():
            all_bin_hits = all_bin_hits.append(bin_hits_df)
        else:
            all_bin_hits = bin_hits_df

# Write out hits folder
all_bin_hits.to_csv(BIN_COUNTS + "all_bin_hits.tsv", sep = '\t', index = False)


###########################
# Generate counts for hits to each bin
###########################
all_bin_counts = all_bin_hits.groupby(["binID", "proteinName"]).size().reset_index(name = "counts")
all_bin_counts.to_csv(BIN_COUNTS + "all_bin_counts.tsv", sep = '\t', index = False)
