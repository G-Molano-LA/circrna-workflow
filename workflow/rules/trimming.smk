#!/usr/bin/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "ALMOST FINISHED" # requires execution to finish it

################################################################################
# Snakefile to trimming RNA-seq reads.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              :
# Last modification : 01-06-2021
################################################################################

# VARIABLES
TRI_READS = expand("{path}/{sample}{ext}", path = config["trimming"]["reads"],
    sample = SAMPLES, ext = [config["trimming"]["suffix"][1], config["trimming"]["suffix"][2]])

TRI_R1 = lambda wildcards: f'{OUTDIR}/trimming/{wildcards.sample}_1_val_1.fq.gz'
TRI_R2 = lambda wildcards: f'{OUTDIR}/trimming/{wildcards.sample}_2_val_2.fq.gz'

R1_TRI      = lambda wildcards: f'{config["trimming"]["reads"]}/{wildcards.sample}{config["trimming"]["suffix"][1]}'
R2_TRI      = lambda wildcards: f'{config["trimming"]["reads"]}/{wildcards.sample}{config["trimming"]["suffix"][2]}'

# TARGET RULE
rule trimming_results:
    input:
        f'{OUTDIR}/quality_control/trimmed_data/summary/multiqc_report.html',
        directory(f'{OUTDIR}/trimming/reports/')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TRIM_GALORE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rule trimming:
    input:
        RAW_READS if QUALITY else TRI_READS
    output:
        reads = expand("{outdir}/trimming/{sample}_{replicate}.fq.gz",  outdir = OUTDIR,
            sample = SAMPLES, replicate = ["1_val_1", "2_val_2"]),
        txt   = expand("{outdir}/trimming/{sample}_{replicate}.fastq.gz_trimming_report.txt",
            outdir = OUTDIR, sample = SAMPLES, replicate = [1,2])
    params:
        outdir = f'{OUTDIR}/trimming/'
    conda: config["envs"]["quality_control"]
    threads: config["trimming"]["threads"]
    priority: 3
    shell:
        """
        # Install Trim Galore
        curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz
        tar xvzf trim_galore.tar.gz --directory=tools/
        rm trim_galore.tar.gz

        # Running
        ./tools/TrimGalore-0.6.6/trim_galore --paired --cores {threads} --output_dir {params.outdir} {input}
        """
rule mv_trimming_reports:
    input:
        txt   = expand("{outdir}/trimming/{sample}_{replicate}.fastq.gz_trimming_report.txt",
            outdir = OUTDIR, sample = SAMPLES, replicate = [1,2])
    output:
        directory(f'{OUTDIR}/trimming/reports/')
    priority: 4
    shell:
        "mv {input} {output}"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FASTQC_POST-TRIMMING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule fastqc2:
    input:
        reads = expand("{path}/trimming/{sample}_{replicate}.fq.gz", path = OUTDIR,
            sample = SAMPLES, replicate = ["1_val_1", "2_val_2"])
    output:
        html = expand("{path}/quality_control/trimmed_data/{sample}_{replicate}_fastqc.html",
            path = OUTDIR, sample = SAMPLES, replicate = ["1_val_1", "2_val_2"]),
        zip  = expand("{path}/quality_control/trimmed_data/{sample}_{replicate}_fastqc.zip",
            path = OUTDIR, sample = SAMPLES, replicate = ["1_val_1", "2_val_2"])
    params:
        outdir = f'{OUTDIR}/quality_control/trimmed_data/'
    conda: config["envs"]["quality_control"]
    # message:
    #     "Starting quality analysis control with FASTQC programm on the "
    #     "following files {input.reads}. The number of threads used are {threads}."
    threads: config["quality_control"]["threads"]
    priority:5
    shell:
        "fastqc -t {threads} {input.reads} --outdir={params.outdir}"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MULTIQC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule multiqc2:
    input:
        zip = expand("{path}/quality_control/trimmed_data/{sample}_{replicate}_fastqc.zip",
            path = OUTDIR, sample = SAMPLES, replicate = ["1_val_1", "2_val_2"])
    output:
        html = f'{OUTDIR}/quality_control/trimmed_data/summary/multiqc_report.html',
    params:
        replace_old = "--force", # revisar que no remplaze al anterior
        outdir      = f'{OUTDIR}/quality_control/trimmed_data/summary/'
    conda: config["envs"]["quality_control"]
    priority: 6
    shell:
        "multiqc --interactive {params.replace_old} {input} --outdir {params.outdir}"
