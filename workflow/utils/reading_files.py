#!bin/bash/python3

###############################################################################
## Python script to obtain filenames.
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Date: 30-06-2021
## Author: G.Molano, LA (gonmola@hotmail.es)
###############################################################################


data= [] # creating an empty list
with open("data/seqs.txt") as names:
    for line in names:
        data.append(line.rstrip('\n'))
SAMPLES=data
