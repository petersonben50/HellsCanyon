#!/bin/sh


cd /Users/benjaminpeterson/Documents/research/HellsCanyon/3D/Brownlee_sampling_sites_v2/original_files/images
ls *jp2 | cut -d "_" -f2 | sort | uniq > image_cluster_list.txt



cd /Users/benjaminpeterson/Documents/research/HellsCanyon/3D/Brownlee_sampling_sites_v2/original_files/
input_file='images/m_4311501_sw_11_1_20150619_20160104.jp2'



####---------------------------------####
# Convert and reduce files
####---------------------------------####
ls images/*jp2 | while read input_file
do
  output_file=`echo $input_file | sed 's/images/images_reduced_size/' | sed 's/jp2/tif/'`
  echo "Reducing" $input_file "to" $output_file
  gdal_translate -of GTiff \
                 -outsize 863 1200 \
                 -r average \
                 $input_file \
                 $output_file
done



####---------------------------------####
# Combine all raster files
####---------------------------------####
threeD_mapping=/Users/benjaminpeterson/Documents/research/HellsCanyon/code/3D_mapping
python $threeD_mapping/merging_image_files.py
