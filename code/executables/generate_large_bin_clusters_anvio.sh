#!/bin/sh

#########################
# executables/generate_large_bin_clusters_anvio.sh
# Benjamin D. Peterson

# This script
#########################


##########################
# Activate anvio environment
##########################
cd $anvioFolder


##########################
# Activate anvio environment
##########################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=''


##########################
# Identify needed number of bins
##########################
anvi-display-contigs-stats $assembly.db \
                            --report-as-text \
                            -o original_summaries/$assembly\_summary.txt
eukBins=`awk -F '\t' -v euk="eukarya (Protista_83)" '$1 == euk { print $2 }' original_summaries/$assembly\_summary.txt`
prokBins=`awk -F '\t' -v prok="bacteria (Bacteria_71)" '$1 == prok { print $2 }' original_summaries/$assembly\_summary.txt`
archBins=`awk -F '\t' -v arch="archaea (Archaea_76)" '$1 == arch { print $2 }' original_summaries/$assembly\_summary.txt`
predictedBins=$(expr $eukBins + $prokBins + $archBins)
requestedBins=$(expr $predictedBins / 3)
echo $requestedBins


##########################
# Run CONCOCT using selected number of bins
##########################
anvi-cluster-contigs -p $assembly.merged/PROFILE.db \
                      -c $assembly.db \
                      -C CONCOCT \
                      -T 16 \
                      --driver concoct \
                      --clusters $requestedBins \
                      --length-threshold 2000 \
                      --just-do-it
