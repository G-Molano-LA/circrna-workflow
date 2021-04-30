#!/usr/bin/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "ALMOST FINISHED" # requires execution to finish it

# Depend on a minimum snakemake version
from snakemake.utils import min_version
min_version("6.0")

# CONFIG FILE
configfile: "src/utils/config.yaml"

# VARIABLES
SAMPLES      = config["samples"]
PATH_raws    = config["path"]["raw_reads"]
PATH_trimmed = config["path"]["trimmed_reads"]
OUTDIR       = config["path"]["outdir"]

# TARGET RULE
rule results:
    input:
        "data/raw_data/samples/multi_report.html",
        "{outdir}/multi_report_trimmed.html"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FASTQC_PRE-TRIMMING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule fastqc1:
    input:
        read = expand("{path}/{sample}{ext}", path = PATH_raws, sample = SAMPLES,
            ext = [config["suffix"]["raw"][1], config["suffix"]["raw"][2]])
    output:
        html = expand("{outdir}/quality_control/raw/{sample}_{replicate}_fastqc.html",
            outdir = OUTDIR, sample = SAMPLES, replicate = [1,2]),
        zip  = expand("{outdir}/quality_control/raw/{sample}_{replicate}_fastqc.zip",
            outdir = OUTDIR, sample = SAMPLES, replicate = [1,2])
    params:
        outdir = f'{OUTDIR}/quality_control/raw/'
    threads: config["threads"]["fastqc"]
    # message:
    #     "Starting quality analysis control with FASTQC programm on the "
    #     "following files {input.reads}. Number of threads used are {threads}."
    priority: 10
    shell:
        "fastqc -t {threads} {input.read2} --outdir={params.outdir}"

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
    priority: 9
    shell:
        "multiqc {params.pdf} {params.replace_old} {input.zip} --outdir {params.outdir}"

#~~~~~~~~~~~~~~~~~~~~TRIM_GALORE_&_FASTQC_POST-TRIMMING~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule trimming:
    input:
        read1 = expand("{path}/{sample}{ext}", path = PATH_raws, sample = SAMPLES,
            ext = config["suffix"]["raw"][1]),
        read2 = expand("{path}/{sample}{ext}", path = PATH_raws, sample = SAMPLES,
            ext = config["suffix"]["raw"][2])
    output:
        reads = expand("{outdir}/{sample}_{replicate}.fq.gz",  outdir = OUTDIR,
            sample = SAMPLES, replicate = ["1_val_1", "2_val_2"]),
        txt   = expand("{outdir}/{sample}_{replicate}.fastq.gz_trimming_report.txt",
            outdir = OUTDIR, sample = SAMPLES, replicate = [1,2])
    params:
        outdir = f'{OUTDIR}/{PATH_trimmed}'
    priority: 8
    shell:
        "trim_galore --paired --output_dir {params.outdir} {input.read1} {input.read2}"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FASTQC_POST-TRIMMING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule fastqc2:
    input:
        reads = expand("{path}/{sample}_{replicate}.fq.gz", path = PATH_trimmed,
            sample = SAMPLES, replicate = ["1_val_1", "2_val_2"])
    output:
        html = expand("{path}/quality_control/trimmed/{sample}_{replicate}_fastqc.html",
            path = OUTDIR, sample = SAMPLES, replicate = ["1_val_1", "2_val_2"]),
        zip  = expand("{path}/quality_control/trimmed/{sample}_{replicate}_fastqc.zip",
            path = OUTDIR, sample = SAMPLES, replicate = ["1_val_1", "2_val_2"])
    params:
        outdir = f'{OUTPUT}/quality_control/trimmed/'
    # message:
    #     "Starting quality analysis control with FASTQC programm on the "
    #     "following files {input.reads}. The number of threads used are {threads}."
    threads: config["threads"]["fastqc"]
    priority: 7
    shell:
        "fastqc -t {threads} {input.reads} --outdir={params.outdir}"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MULTIQC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule multiqc2:
    input:
        zip = expand("{path}/quality_control/trimmed/{sample}_{replicate}_fastqc.zip",
            path = OUTDIR, sample = SAMPLES, replicate = ["1_val_1", "2_val_2"])
    output:
        html = f'{OUTDIR}/quality_control/trimmed/multi_report_trimmed.html',
        pdf  = f'{OUTDIR}/quality_control/trimmed/multi_report_trimmed.pdf'
    params:
        pdf         = "--pdf",
        replace_old = "--force", # revisar que no remplaze al anterior
        outdir      = f'{OUTDIR}/quality_control/trimmed/'
    priority: 6
    shell:
        "multiqc {params.pdf} {params.replace_old} {input} --outdir {params.outdir}"
