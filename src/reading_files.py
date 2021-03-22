#!bin/bash/python3

###############################################################################
## This file is an attempt to obtain filenames
###############################################################################

data= [] # creating an empty list
with open("files.txt") as fp:
    for line in fp:
        data.append(line.strip())
print(data)
