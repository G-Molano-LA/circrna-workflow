#!/bin/bash


###############################################################################
## Obtainning filenames of data samples
## Author: G. Molano, LA
###############################################################################
cd data/trimmed_data/ # data directory

# creating files var with all output of ls
files=( $(ls * ))
declare -a arr

# passing command-line argument
ext=$1

# Checking if argument has been supplied
if [ -z "$ext" ]; then
  echo "ERROR: No argument supplied. Please specify an expression to delete."
else
  # Deleting expression from filenames
  for name in "${files[@]}"
    do
      arr+=(${name%%$ext})
    done
  # Deleting duplicates generated from forward and reverse samples
  filenames=($(printf "%s\n" "${arr[@]}" | sort -u))
  echo "The generated samaple names are stored in 'seqs.txt' file. Sample names: ${filenames[@]}"

fi

cd .. # to data/ directory
# Redirecting output + Deleting duplicates generated from forward and reverse samples
printf "%s\n" "${filenames[@]}"> seqs.txt

cd .. # to circrna_workflow/ directory
