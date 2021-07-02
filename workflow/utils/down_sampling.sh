#!/bin/bash

###############################################################################
## Script to down-sampling trimmed data
## Author: G. Molano, LA (gonmola@hotmail.es)
###############################################################################

# 1. Installation of down-sampling tool for processing sequences in FASTQ formats
#    (source: https://github.com/lh3/seqtk)
  # git clone https://github.com/lh3/seqtk.git;
  # cd seqtk; make

cd data/trimmed_data # data directory
files=( $(ls * )) # creating files var with all output of ls

for filename in "${files[@]}" # accessing all content of files var
do
  newname="${filename/_/_sub_}" # substitution of the first underscore matched
  echo "Down-sampling $filename ..."
  seqtk sample -s100 $filename 10000 > $newname # Down-sampling
  echo "Completed."
done

cd ..; cd ..
