#####################################
# code/generalUse/calculate_depth_contigs_in_bins.py
# Benjamin D. Peterson
#####################################

#####################################
# Load up libraries
#####################################
import os
import sys
import pandas as pd

#####################################
# Read in input from command line
#####################################
depthTableFile = sys.argv[1]
filteredLength = int(sys.argv[2])
scaffolds2binFile = sys.argv[3]
outputFileName = sys.argv[4]


#####################################
# Open up depth file
#####################################
depthTable = pd.read_table(depthTableFile, names = ['contigID', 'locus', 'depth'])

#####################################
# Filter out residues at start of contig
#####################################
depthTable = depthTable[depthTable['locus'] >= filteredLength]


#####################################
# Filter out residues at end of contig
#####################################

# First, find the max contig length
lengthOfContig = depthTable[['contigID', 'locus']].groupby('contigID').max()
lengthOfContig.rename(columns = {'locus':'lengthOfContig'}, inplace = True)
# Save out original length
lengthOfContigOriginal = lengthOfContig
lengthOfContigOriginal
lengthOfContig['maxLengthToInclude'] = lengthOfContig['lengthOfContig'] - filteredLength
# Then, join max contig length DF with depth table
depthTable = pd.merge(depthTable, lengthOfContig, on='contigID', how='outer')
# Filter out end of contig
depthTable = depthTable[depthTable['locus'] <= depthTable['maxLengthToInclude']]
# Keep the depth and contigID
depthTable = depthTable[['contigID', 'depth']]


#####################################
# Sum depth by contig
#####################################
# Find average coverage across contig
depthTableMeanDepth = depthTable.groupby('contigID').mean()


#####################################
# Add in length information
#####################################
depthTableDepthLength = pd.merge(depthTableMeanDepth, lengthOfContig, on='contigID', how='outer')


#####################################
# Read scaffolds2bin file
#####################################
scaffolds2bin = pd.read_table(scaffolds2binFile, names = ['contigID', 'binID'])


#####################################
# Combine the S2B and depth files, then aggregate by bin
#####################################
print("Merging" + depthTableFile + "file with" + scaffolds2binFile)
depthTableBinAll = pd.merge(depthTableDepthLength, scaffolds2bin, on = 'contigID')
depthTableBinAll = depthTableBinAll[['binID', 'lengthOfContig', 'depth']]
depthTableBinAll[['totalCoverage']] = depthTableBinAll['lengthOfContig']*depthTableBinAll['depth']


#####################################
# Aggregate by bin
#####################################
print("Aggregating depth by contig")
depthBinSum = depthTableBinAll.groupby('binID').sum()
depthBinSum['meanCoverageBin'] = depthBinSum['totalCoverage'] / depthBinSum['lengthOfContig']
depthBinFinal = depthBinSum[['meanCoverageBin']]


#####################################
# Save it out
#####################################
depthBinFinal.to_csv(outputFileName)
