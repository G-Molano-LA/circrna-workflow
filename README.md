

### Table of Contents
* [Requirements](#requirements)
  1. [Snakemake: workflow system manager](#snakemake-workflow-system-manager)
  2. [Setting your working directory](#setting-your-working-directory)

***
# Requirements
To execute the CircWokflow, the following programs are required:

* [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) >= 6.1.1
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 0.11.9
* [MultiQC](https://multiqc.info/) 1.10.1
* [TrimGalore](https://github.com/FelixKrueger/TrimGalore) 0.6.4
* [Burrow-Wheeler Aligner](https://github.com/lh3/bwa) 0.7
* [CIRI2](https://sourceforge.net/projects/ciri/files/CIRI2/) 2.0.6
* [CircExplorer2](https://github.com/YangLab/CIRCexplorer2) 2.3.8
* [CIRIquant](https://github.com/bioinfo-biols/CIRIquant) 1.1.2

However, we do not recommended to install any of them manually, whilst handle them automatically by using Conda environments is the best option.

## Snakemake: workflow system manager
Snakemake is the workflow management system used to implement CircWorkflow. It uses a domain specific language based on Python to create transparent workflows, ensuring the reproducibility and automation of data analysis.

### Installation of Dependencies
#### Python3
```{bash}
sudo apt-get update
sudo apt-get install python3.8 python3-pip
```
#### Anaconda
To facilitate the installation and managing of the rest of programmes , we will use the [Anaconda](https://www.anaconda.com/) distribution and the package manager [conda](https://conda.io/projects/conda/en/latest/index.html).

Anaconda is an optimized Python and R distribution, having pre-built and pre-configured collection of packages that can be installed and used on a system. Anaconda uses a package manager named `conda`, which can not only built and manage software from Python language, but also from any type of programming language. Additionally, we also will use another package manager called `mamba` to make easier the installation of `snakemake` tool.

package manager
: is a tool that automates the process of installing, updating, and removing packages.

To install Anaconda distribution on Linux follow the steps. For other operating system check the [man page](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html):
1. Download the installer of Anaconda for Linux ([here](https://www.anaconda.com/products/individual)].
2. Open a terminal window (ALT+Ctl+T) and go to the directory containing the downloaded installer.
3. Run:
```{bash}
bash Anaconda-latest-Linux-x86_64.sh
conda update
```
4. Follow the prompts on the installer screens. In front of any doubt, accept the defaults settings. They can be changed later.
5. To make the changes take effect, close and re-open your terminal window.
6. Test your installation. In your terminal window run the command `conda list`, which will show a list of installed packages if it has been installed correctly.

In order to make visible to your system and be able to run all the pre-installed programmes with Anaconda distribution from anywhere, is necessary to add the path of Anaconda to your `$PATH` bash variable:

1. Open the following document with root privileges:
```{bash}
sudo nano /etc/profile
```
2. Modify the document, adding at the top of the document:
```{bash}
export PATH=$PATH:/your/path/to/anaconda3/bin:/your/path/to/anaconda3/condabin
```
This modification will affect global settings, which means that it will be run upon login for all current and future users. To modify the `$PATH` variable only for a specific user:
```{bash}
sudo nano ~/.profile

# Add at the top of the profile file:
export PATH=$PATH:/your/path/to/anaconda3/bin
```
3. Once the profile file has been modified, reload it and ensure that `$PATH` has added the path of Anaconda distribution.
```{bash}
source /etc/profile
$PATH
```
4. Add repositories (called `channels`) to install packages. The order of adding these channels is important.
```{bash}
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

#### Installing Snakemake
To simplify the installation of this programme, we will use the `mamba` package manager instead of the default `conda` manager.

```{bash}
conda update
conda install mamba
mamba install -c conda-forge -c bioconda snakemake
```

## Setting your working directory

1. Download and unzip this repository.
2. Modify the config file located at `workflow/config/config.yaml`. Specify the workflow module that you want to execute and complete the general and specific configurations.

3. Perform a dry run to ensure that there is no hidden problems.
```
snakemake -n
```

4. After ensuring this, you can now execute the process:
```
snakemake --cores NUMBER --use-conda
```
