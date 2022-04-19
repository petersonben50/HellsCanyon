#!/bin/sh

# mise_en_place.sh
# Benjamin D. Peterson

# This is the script that pulls down the scripts
# from Github that are needed to run this workflow.

cd ~/HellsCanyon/code
git clone https://github.com/petersonben50/homemade_bioinformatic_scripts.git
cd ~/HellsCanyon/code


# Get updates
cd ~/HellsCanyon/code/homemade_bioinformatic_scripts
#git fetch
#git diff ...origin
#git pull
