
__author__ = "G. Molano, LA (gonmola@hotmail.es)"

import subprocess
import yaml
import sys

# Passing arguments
GENOME      = str(sys.argv[1])
PATH_genome = str(sys.argv[2])
OUTFILE     = str(sys.argv[3])

# Capturing conda path
conda_path = subprocess.check_output('which bwa', shell = True)
conda_path = conda_path.rstrip('/bwa\n')

dir_path = subprocess.check_output('pwd', shell = True)
dir_path = dir_path.rstrip('\n')

# Creating dictionary for yaml file
dic_yaml = {'name'      : GENOME,
            'tools'     :
                {'bwa'      : '{}/bwa'.format(conda_path),
                 'hisat2'   : '{}/hisat2'.format(conda_path),
                 'stringtie': '{}/stringtie'.format(conda_path),
                 'samtools' : '{}/samtools'.format(conda_path)},
            'reference' :
                {'fasta'         : '{}/{}/{}.fna'.format(dir_path, PATH_genome, GENOME),
                 'gtf'           : '{}/{}/{}_ann.gtf'.format(dir_path, PATH_genome, GENOME),
                 'bwa_index'     : '{}/{}/bwa/{}.fna'.format(dir_path, PATH_genome, GENOME),
                 'hisat_index' : '{}/{}/hisat2/{}'.format(dir_path, PATH_genome, GENOME)
                 }
            }

# Writing YAML file
with open(OUTFILE, 'w') as file:
    documents = yaml.dump(dic_yaml, file)
print("Done")
