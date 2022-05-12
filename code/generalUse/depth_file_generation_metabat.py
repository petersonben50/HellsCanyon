
# Must have bioinformatics conda environment
# active, or any environment with MetaBat2
# installed.

# BAM files must be in form:
# $metagenome\_to_$assembly.bam

import sys
import os
import subprocess

list_of_metagenomes = sys.argv[1]
mapping_file_location = sys.argv[2]
assembly_name = sys.argv[3]
output_location = sys.argv[4]

#list_of_metagenomes = "/home/GLBRCORG/bpeterson26/HellsCanyon/metadata/lists/2017_analysis_metagenomes.txt"
#mapping_file_location = "/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/2017_analysis_bins/binning/mapping"
#assembly_name = "fall2017cluster1"
#output_location = "/home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/2017_analysis_bins/binning/autoBinning/metabat2"
# Read in list of metagenomes

jgi_command = "jgi_summarize_bam_contig_depths --outputDepth " + output_location + "/depth_to_" + assembly_name + ".txt "
with open(list_of_metagenomes, "r") as bam_list:
  for metagenome in bam_list:
    stripped_line = metagenome.strip()
    jgi_command = jgi_command + mapping_file_location + "/" + stripped_line + "_to_" + assembly_name + ".bam "

print(jgi_command)
subprocess.run(jgi_command, shell=True)
