#!/bin/sh

####---------------------------------####
# code/3D_mapping/merging_image_files.sh
# Benjamin D. Peterson

# This set of scripts reduces the size of
# the image files we retrieved for the 3D
# map of Brownlee and then merges them
# using gdal_merge.
####---------------------------------####


####---------------------------------####
# Get set up
####---------------------------------####
conda activate py_viz
cd /Users/benjaminpeterson/Documents/research/HellsCanyon/3D/Brownlee_sampling_sites_v2/original_files/



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
gdalwarp -t_srs EPSG:4269 images_merged.tif images_merged_4269.tif
