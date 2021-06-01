#!/usr/bin/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "ALMOST FINISHED" # requires execution to finish it
################################################################################
# Snakefile to realize a quality control of RNA-seq reads.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              :
# Last modification : 26-05-2021
################################################################################
RAW_READS = expand("{path}/{sample}{ext}", path = config["quality_control"]["reads"], sample = SAMPLES,
    ext = [config["quality_control"]["suffix"][1],config["quality_control"]["suffix"][2]] )

# TARGET RULE
rule quality_control_results:
    input:
        html = f'{OUTDIR}/quality_control/raw_data/summary/multiqc_report.html'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FASTQC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule fastqc1:
    input:
        RAW_READS
    output:
        html = expand("{outdir}/quality_control/raw_data/{sample}_{replicate}_fastqc.html",
            outdir = OUTDIR, sample = SAMPLES, replicate = [1,2]),
        zip  = expand("{outdir}/quality_control/raw_data/{sample}_{replicate}_fastqc.zip",
            outdir = OUTDIR, sample = SAMPLES, replicate = [1,2])
    params:
        outdir = f'{OUTDIR}/quality_control/raw_data/'
    threads: config["trimming"]["threads"]
    conda: config["envs"]["quality_control"]
    # message:
    #     "Starting quality analysis control with FASTQC programm on the "
    #     "following files {input.reads}. Number of threads used are {threads}."
    priority: 10
    shell:
        "fastqc -t {threads} {input} --outdir={params.outdir}"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MULTIQC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule multiqc1:
    input:
        zip  = expand("{outdir}/quality_control/raw_data/{sample}_{replicate}_fastqc.zip",
            outdir = OUTDIR, sample = SAMPLES, replicate = [1,2])
    output:
        html = f'{OUTDIR}/quality_control/raw_data/summary/multiqc_report.html',
    params:
        replace_old = "--force", # revisar que no remplaze al anterior
        outdir      = f'{OUTDIR}/quality_control/raw_data/summary/'
    conda: config["envs"]["quality_control"]
    priority: 9
    shell:
        "multiqc --interactive {params.replace_old} {input.zip} --outdir {params.outdir}"
