#!bin/bash/python3

###############################################################################
## This file is an attempt to obtain filenames
###############################################################################


data= [] # creating an empty list
with open("data/seqs.txt") as fp:
    for line in fp:
        data.append(line.rstrip('_[1,2]_val_[1,2].fq.gz\n'))
data=list(dict.fromkeys(data)) # removing duplicates
print(data)
