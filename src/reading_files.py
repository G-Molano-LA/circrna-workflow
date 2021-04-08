#!bin/bash/python3

###############################################################################
## This file is an attempt to obtain filenames
###############################################################################


data= [] # creating an empty list
with open("data/seqs.txt") as names:
    for line in names:
        data.append(line.rstrip('\n'))
SAMPLES=data
