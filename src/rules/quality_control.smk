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
RAW_READ1 = expand("{path}/{sample}{ext}", path = config["quality_control"]["reads"], sample = SAMPLES,
    ext = config["quality_control"]["suffix"][1])
RAW_READ2 = expand("{path}/{sample}{ext}", path = config["quality_control"]["reads"], sample = SAMPLES,
    ext = config["quality_control"]["suffix"][2])

# TARGET RULE
rule quality_control_results:
    input:
        f'{OUTDIR}/quality_control/raw/multi_report.html'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FASTQC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule fastqc1:
    input:
        read1 = RAW_READ1,
        read2 = RAW_READ2
    output:
        html = expand("{outdir}/quality_control/raw/{sample}_{replicate}_fastqc.html",
            outdir = OUTDIR, sample = SAMPLES, replicate = [1,2]),
        zip  = expand("{outdir}/quality_control/raw/{sample}_{replicate}_fastqc.zip",
            outdir = OUTDIR, sample = SAMPLES, replicate = [1,2])
    params:
        outdir = f'{OUTDIR}/quality_control/raw/'
    threads: config["trimming"]["threads"]
    conda: config["envs"]["quality_control"]
    # message:
    #     "Starting quality analysis control with FASTQC programm on the "
    #     "following files {input.reads}. Number of threads used are {threads}."
    priority: 10
    shell:
        "fastqc -t {threads} {input.read1} {input.read2} --outdir={params.outdir}"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MULTIQC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule multiqc1:
    input:
        zip  = expand("{outdir}/quality_control/raw/{sample}_{replicate}_fastqc.zip",
            outdir = OUTDIR, sample = SAMPLES, replicate = [1,2])
    output:
        html = f'{OUTDIR}/quality_control/raw/multi_report.html',
        pdf  = f'{OUTDIR}/quality_control/raw/multi_report.pdf'
    params:
        pdf         = "--pdf",
        replace_old = "--force", # revisar que no remplaze al anterior
        outdir      = f'{OUTDIR}/quality_control/raw/'
    conda: config["envs"]["quality_control"]
    priority: 9
    shell:
        "multiqc {params.pdf} {params.replace_old} {input.zip} --outdir {params.outdir}"
