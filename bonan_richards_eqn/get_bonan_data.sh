#!/bin/bash

# This script downloads a specific version of a Matlab script from a GitHub repository
# and generates data files for two different soil types using that script.
# The output data files are saved in the paths bonanmodeling/sp_08_01/bonan_data_clay.txt
# and bonanmodeling/sp_08_01/bonan_data_sand.txt, relative to the directory
# where this script is run.

# Parse the input argument for the output directory
output_dir=$1

# Clone the repository containing the Matlab script we need
git clone https://github.com/gbonan/bonanmodeling.git
cd bonanmodeling

# Checkout the commit that we know has what we need
git -c advice.detachedHEAD=false checkout a10cf764013be58c2def1dbe7c7e52a3213e061e

# Load the Matlab module on Caltech's Central cluster
# NOTE: You will need to change it if running on a different machine
module load matlab/r2024a

# Run the Matlab script to generate the data for the clay soil type
matlab -batch "run('sp_08_01/sp_08_01.m'); exit;"

# Remove the header (first line) and rename the generated data file
tail -n +2 sp_08_01/data1.txt > ../$output_dir/bonan_data_clay.txt

# Comment out lines 60-66 of sp_08_01.m and uncomment lines 69-75
# This changes the soil type parameters used for the experiment
awk 'NR>59&&NR<67{$1 = "%" substr($0, 2)}' sp_08_01/sp_08_01.m
awk 'NR>68&&NR<76{$1 = " " substr($0, 2)}' sp_08_01/sp_08_01.m

# Run the Matlab script again to generate the data for the sand soil type
matlab -batch "run('sp_08_01/sp_08_01.m'); exit;"

# Remove the header (first line) and rename the generated data file
tail -n +2 sp_08_01/data1.txt > ../$output_dir/bonan_data_sand.txt
