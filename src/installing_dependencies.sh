#!/bin/bash

################################################################################
## Dependencies Installation on Linux machines
################################################################################

# 1.Setting your working directory
mkdir -p circrna_workflow/raw_data circrna_workflow/ciri circrna_workflow/src
cd circrna_workflow/


# 2. Installing CIRI2
wget https://sourceforge.net/projects/ciri/files/CIRI2/CIRI_v2.0.6.zip
unzip CIRI_v2.0.6.zip
rm CIRI_v2.0.6.zip
mv CIRI_v2.0.6/CIRI2.pl ciri/
rm CIRI2_v2.0.6

# 3. Installing Anaconda
wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
bash Anaconda3-2020.11-Linux-x86_64.sh
rm Anaconda3-2020.11-Linux-x86_64.sh
bash --login

## 3.1. Conda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# 4. Installing snakemake
conda install mamba
mamba install -c conda-forge -c bioconda snakemake

# 5. Creating a conda environment
conda env list
wget https://github.com/G-Molano-LA/circrna-workflow/raw/main/configuration.yml
conda create --name circrna_env configuration.yml
conda env list
conda activate circrna_env
