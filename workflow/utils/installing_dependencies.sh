#!/bin/bash

################################################################################
## Dependencies Installation on Linux machines
################################################################################

# 1.Setting your working directory
mkdir -p circrna_workflow/data/raw_data7samples circrna_workflow/src \
 circrna_workflow/data/mapped_data
cd circrna_workflow/

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
conda update
conda install mamba -c conda-forge
mamba install -c conda-forge -c bioconda snakemake
