
__author__ = "G. Molano, LA (gonmola@hotmail.es)"

import subprocess
import yaml
import sys

# Passing arguments
GENOME      = str(sys.argv[1])
PATH_genome = str(sys.argv[2])
OUTFILE     = f'{str(sys.argv[3])}/ciriquant_config.yaml'

# Capturing shell output
conda_path = subprocess.run('snakemake --use-conda --list-conda-envs',
                        shell = True, text = True, capture_output = True)
start_index = conda_path.stdout.find('snakemake')
end_index = len(conda_path.stdout)
conda_path = conda_path.stdout[start_index:end_index]
conda_path = conda_path.rstrip('\n')

dir_path = subprocess.run('pwd', shell = True, text = True, capture_output = True)
dir_path = dir_path.stdout.rstrip('\n')

# Creating dictionary for yaml file
dic_yaml = {'name'      : GENOME,
            'tools'     :
                {'bwa'      : f'{dir_path}/{conda_path}/bin/bwa',
                 'hisat2'   : f'{dir_path}/{conda_path}/bin/hisat2',
                 'stringtie': f'{dir_path}/{conda_path}/bin/stringtie',
                 'samtools' : f'{dir_path}/{conda_path}/bin/samtools'},
            'reference' :
                {'fasta'         : f'{PATH_genome}/{GENOME}.fna',
                 'gft'           : f'{PATH_genome}/{GENOME}_ann.gtf',
                 'bwa_index'     : f'{PATH_genome}/bwa/{GENOME}',
                 'hitsat2_index' : f'{PATH_genome}/hisat2/{GENOME}'
                 }
            }

# Writing YAML file
with open(OUTFILE, 'w') as file:
    documents = yaml.dump(dic_yaml, file)
print("Done")
