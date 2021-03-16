#!/usr/bin/python3

# Depend on a minimum snakemake version
from snakemake.utils import min_version
min_version("6.0")


## 1. QUALITY CONTROL #########################################################
rule quality_control:
    input:
        reads=expand("raw_data/{sample}_{replicate}.{ext}", replicate=[1,2],
        ext=["fastq.gz", "fastq"])
    output:
        html="fastqc/{sample}.html"
        zip="fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary
                                        #+ multiqc find the file.
    log:
        "fastqc/logs/{sample}_fastqc.log"
    threads: 3
    message:
    "Starting quality analysis control with FASTQC programm on the "
    "following files {input.reads}. Number of threads used are {threads}."
    shell:
        "fastqc -t {threads} {input.reads} 2>{log}"

rule multiqc_report:
    input:
        "fastqc/{sample}_fastqc.zip"
    output:
        html="raw_data/fastqc/multi_report.html"
        pdf="raw_data/fastqc/multi_report.pdf"
    params:
        pdf="--pdf"
        replace_old="--force"
    log:
        "fastqc/logs/{sample}_multiqc.log"
    shell:
        "multiqc {params.pdf} {params.replace_old} {input} 2>{log}"

rule trimming:
    input:
        reads=expand("raw_data/{sample}_{replicate}.{ext}", replicate=[1,2],
        ext=["fastq.gz", "fastq"])
    output:
        reads=expand("trimmed_data/{sample}.{replicate}.fq.gz",
            replicate=["1_val_1", "2_val_2"])
        txt=expand("trimmed_data/{sample}.{replicate}.fastq.gz_trimming_report.txt",
            replicate=[1,2])
    shell:
        "trim_galore --paired -o trimed_data {input}"

rule quality_control_2:
    input:
        reads=expand("trimmed_data/{sample}.{replicate}.fq.gz",
            replicate=["1_val_1", "2_val_2"])
    output:
        html="trimmed_data/{sample}.html"
        zip="trimmed_data/{sample}_fastqc.zip"
    log:
        "trimmed_data/logs/{sample}_fastqc.log"
    threads: 3
    message:
    "Starting quality analysis control with FASTQC programm on the "
    "following files {input.reads}. The number of threads used are {threads}."
    shell:
        "fastqc -t {threads} {input.reads} 2>{log}"

rule multiqc_report_2:
    input:
        "trimmed_data/{sample}_fastqc.zip"
    output:
        html="trimmed_data/multi_report.html"
        pdf="trimmed_data/multi_report.pdf"
    params:
        pdf="--pdf"
        replace_old="--force"
    log:
        "trimmed_data/logs/{sample}_multiqc.log"
    shell:
        "multiqc {params.pdf} {params.replace_old} {input} 2>{log}"
