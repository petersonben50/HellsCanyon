
# Protocol built from here:
# https://www.neonscience.org/resources/learning-hub/tutorials/merge-lidar-geotiff-py

####---------------------------------####
# Load up needed libraries
####---------------------------------####

import numpy as np
import matplotlib.pyplot as plt
import os, glob
from osgeo import gdal


####---------------------------------####
# Set up file names
####---------------------------------####
# Set up an argument parser
"""
parser = argparse.ArgumentParser()
parser.add_argument('--file_names')
parser.add_argument('--input_image_files')
parser.add_argument('--output_image_files')

inputs = parser.parse_args()
FILE_NAMES = inputs.file_names
INPUT_IMAGE_LOCATION = inputs.input_image_location
OUTPUT_IMAGE_FILES = inputs.output_image_files
"""

FILE_NAMES = ''
INPUT_IMAGE_LOCATION = '/Users/benjaminpeterson/Documents/research/HellsCanyon/3D/Brownlee_sampling_sites_v2/original_files/images_reduced_size'
OUTPUT_IMAGE_LOCATION = '/Users/benjaminpeterson/Documents/research/HellsCanyon/3D/Brownlee_sampling_sites_v2/merged_files'
OUTPUT_FILE_NAME = "images_merged.tif"
OUTPUT_IMAGE_FILE = OUTPUT_IMAGE_LOCATION + "/" + OUTPUT_FILE_NAME


####---------------------------------####
# Read in variables
####---------------------------------####
files_to_mosaic = glob.glob(INPUT_IMAGE_LOCATION + "/*.tif")
files_string = " ".join(files_to_mosaic)

command = "gdal_merge.py -co 'ALPHA=NO' -o " + OUTPUT_IMAGE_FILE + " " + files_string
print(os.popen(command).read())
