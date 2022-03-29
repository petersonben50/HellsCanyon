#!/bin/sh

######################
# code/binning/binAnalysis_depth.sh
# Benjamin D. Peterson

# This set of scripts contains the
# scripts needed to pull out the depth
# data for the
######################


####################################################
####################################################
# Pull out coverage data for hgcA+ bins
####################################################
####################################################

binAnalysis=~/HellsCanyon/dataEdited/binAnalysis
mkdir $binAnalysis/depth
cd /home/GLBRCORG/bpeterson26/HellsCanyon/code/
chmod +x executables/aggregate_depth_bins.sh
condor_submit submission/aggregate_depth_bins.sub



####################################################
####################################################
# Pull out coverage data for HQ bins
####################################################
####################################################
cd /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binning/autoBinning/hqBinSet/
mkdir depth
cd /home/GLBRCORG/bpeterson26/HellsCanyon/code/
condor_submit submission/aggregate_depth_bins_hqBins.sub
