#!bin/bash/python3

###############################################################################
## This file is an attempt to obtain filenames
###############################################################################


data= [] # creating an empty list
with open("files.txt") as fp:
    for line in fp:
        data.append(line.rstrip('_[1,2].fastq.gz\n'))
data=list(dict.fromkeys(data))
print(data)
