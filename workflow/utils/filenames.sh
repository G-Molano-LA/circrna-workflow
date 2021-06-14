#!/bin/bash


###############################################################################
## Obtainning filenames of data samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              :
# Last modification : 01-06-2021
################################################################################
# [USAGE]:
# ./utils/filenames.sh <data_directory> <filenames_extension>


# passing command-line argument
dir=$1
ext=$2

pushd $dir

# creating files var with all output of ls
files=( $(ls ))
declare -a arr

# Checking if argument has been supplied
if [ -z "$ext" ]; then
  echo "ERROR: No argument supplied. Please specify an expression to delete."
else
  # Deleting expression from filenames
  for name in "${files[@]}"
    do
      arr+=("${name%%$ext}")
    done
  # Deleting duplicates generated from forward and reverse samples
  filenames=($(printf "%s\n" "${arr[@]}" | sort -u))
  echo "The generated samaple names are stored in 'data/seqs.txt' file. Sample names: ${filenames[@]}"

fi

popd

# Redirecting output + Deleting duplicates generated from forward and reverse samples
printf "%b\n- " "- ${filenames[@]}" > seqs.txt
