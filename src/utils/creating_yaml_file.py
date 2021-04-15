
__author__ = "G. Molano, LA (gonmola@hotmail.es)"

import subprocess
import yaml

# Capturing shell output
conda_path = subprocess.run('snakemake --snakefile src/main/quantification.Snakefile \
                                   --use-conda --list-conda-envs',
                        shell = True, text = True, capture_output = True)
start_index = conda_path.stdout.find('snakemake')
end_index = len(conda_path.stdout)
conda_path = conda_path.stdout[start_index:end_index]
conda_path = conda_path.rstrip('\n')

dir_path = subprocess.run('pwd', shell = True, text = True, capture_output = True)
dir_path = dir_path.stdout.rstrip('\n')

# Creating dictionary for yaml file
dic_yaml = {'name'      : 'GRCh38',
            'tools'     :
                {'bwa'      : f'{dir_path}/{conda_path}{/bin/bwa}',
                 'hisat2'   : f'{dir_path}/{conda_path}{/bin/hisat2}',
                 'stringtie': f'{dir_path}/{conda_path}{/bin/stringtie}',
                 'samtools' : f'{dir_path}/{conda_path}{/bin/samtools}'},
            'reference' :
                {'fasta'         : f'{dir_path}/data/raw_data/GRCh38.fa',
                 'gft'           : f'{dir_path}/data/raw_data/GRCh38_ann.gtf',
                 'bwa_index'     : f'{dir_path}/data/raw_data/bwa/GRCh38',
                 'hitsat2_index' : f'{dir_path}/data/raw_data/hisat2/GRCh38'}
                 }

# Writing YAML file
with open(r'libs/ciriquant/config_file.yaml', 'w') as file:
    documents = yaml.dump(dic_yaml, file)
print("Done")
