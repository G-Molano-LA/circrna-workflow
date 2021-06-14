#!/bin/python2

__author__ = "G. Molano, LA (gonmola@hotmail.es)"

import subprocess
import yaml
import sys

# Passing arguments
GENOME       = str(sys.argv[1])
HISAT2_INDEX = str(sys.argv[2])
BWA_INDEX    = str(sys.argv[3])
REFERENCE    = str(sys.argv[4])
ANNOTATION   = str(sys.argv[5])
OUTFILE      = str(sys.argv[6])

HISAT2_BOL     = str(sys.argv[7])
BWA_BOL        = str(sys.argv[8])
REFERENCE_BOL  = str(sys.argv[9])
ANNOTATION_BOL = str(sys.argv[10])

# Capturing conda path
conda_path = subprocess.check_output('which bwa', shell = True)
conda_path = conda_path.rstrip('/bwa\n')

dir_path = subprocess.check_output('pwd', shell = True)
dir_path = dir_path.rstrip('\n')

# Conditional path assignment
paths = [ REFERENCE, ANNOTATION, BWA_INDEX, HISAT2_INDEX ]
files = [ REFERENCE_BOL, ANNOTATION_BOL, BWA_BOL, HISAT2_BOL ]

count = -1
for boolean in files:
    count = count + 1
    if boolean is True:
        pass
    else:
        paths[count] = '{}/{}'.format(dir_path, paths[count])
# Creating dictionary for yaml file
dic_yaml = {'name'      : GENOME,
            'tools'     :
                {'bwa'      : '{}/bwa'.format(conda_path),
                 'hisat2'   : '{}/hisat2'.format(conda_path),
                 'stringtie': '{}/stringtie'.format(conda_path),
                 'samtools' : '{}/samtools'.format(conda_path)},
            'reference' :
                {'fasta'         : paths[0],
                 'gtf'           : paths[1],
                 'bwa_index'     : paths[2],
                 'hisat_index'   : paths[3]
                 }
            }

# Writing YAML file
with open(OUTFILE, 'w') as file:
    documents = yaml.dump(dic_yaml, file)
print("Done")
