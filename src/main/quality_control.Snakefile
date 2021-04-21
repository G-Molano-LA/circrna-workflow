#!/usr/bin/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "ALMOST FINISHED" # requires execution to finish it

# Depend on a minimum snakemake version
from snakemake.utils import min_version
min_version("6.0")

###############################################################################
data= [] # creating an empty list
with open("data/seqs.txt") as fp:
    for line in fp:
        data.append(line.rstrip('_[1,2].fastq.gz\n'))
data=list(dict.fromkeys(data)) # removing duplicates
SAMPLES=data
###############################################################################
rule results:
    input:
        "data/raw_data/samples/multi_report.html",
        "data/trimmed_data/multi_report_trimmed.html"


"""
You need to have at least one rule (target rule) that does not produce any output,
but takes as input all your expected output. Snakemake will take that rule as
first rule, and then checks which rules produce the input of the rule. Then it
checks, which rules produce the inputs for those rules, etc.
"""

## 1. QUALITY CONTROL #########################################################
rule quality_control:
    input:
        reads=expand("data/raw_data/samples/{sample}_{replicate}.fastq.gz",
            sample=SAMPLES, replicate=[1,2])
    output:
        html=expand("data/raw_data/samples/{sample}_{replicate}_fastqc.html",
            sample=SAMPLES, replicate=[1,2]),
        zip=expand("data/raw_data/samples/{sample}_{replicate}_fastqc.zip",
            sample=SAMPLES, replicate=[1,2])
    log:
        "logs/fastqc/fastqc_1.log"
    threads: 3
    # message:
    #     "Starting quality analysis control with FASTQC programm on the "
    #     "following files {input.reads}. Number of threads used are {threads}."
    priority: 10
    shell:
        "fastqc -t {threads} {input.reads} 2>{log}"

rule multiqc_report:
    input:
        expand("data/raw_data/samples/{sample}_{replicate}_fastqc.zip",
                sample=SAMPLES, replicate=[1,2])
    output:
        html="data/raw_data/samples/multi_report.html",
        pdf="data/raw_data/samples/multi_report.pdf"
    params:
        pdf="--pdf",
        replace_old="--force" # revisar que no remplaze al anterior
    log:
        "logs/multiqc/multiqc_1.log"
    priority: 9
    shell:
        "multiqc {params.pdf} {params.replace_old} {input} 2>{log}"

rule trimming:
    input:
        reads=expand("data/raw_data/samples/{sample}_{replicate}.fastq.gz",
            sample=SAMPLES, replicate=[1,2])
    output:
        reads=expand("data/trimmed_data/{sample}_{replicate}.fq.gz",
             sample=SAMPLES, replicate=["1_val_1", "2_val_2"]),
        txt=expand("data/trimmed_data/{sample}_{replicate}.fastq.gz_trimming_report.txt",
             sample=SAMPLES, replicate=[1,2])
    priority: 8
    shell:
        "trim_galore --paired -o trimed_data {input}"

rule quality_control_2:
    input:
        reads=expand("data/trimmed_data/{sample}_{replicate}.fq.gz",
             sample=SAMPLES, replicate=["1_val_1", "2_val_2"])
    output:
        html=expand("data/trimmed_data/{sample}_{replicate}_fastqc.html",
             sample=SAMPLES, replicate=["1_val_1", "2_val_2"]),
        zip=expand("data/trimmed_data/{sample}_{replicate}_fastqc.zip",
             sample=SAMPLES, replicate=["1_val_1", "2_val_2"])
    log:
        "logs/fastqc/fastqc_2.log"
    # message:
    #     "Starting quality analysis control with FASTQC programm on the "
    #     "following files {input.reads}. The number of threads used are {threads}."
    priority: 7
    shell:
        "fastqc -t {threads} {input.reads} 2>{log}"

rule multiqc_report_2:
    input:
        zip=expand("data/trimmed_data/{sample}_{replicate}_fastqc.zip",
             sample=SAMPLES, replicate=["1_val_1", "2_val_2"])
    output:
        html="data/trimmed_data/multi_report_trimmed.html",
        pdf="data/trimmed_data/multi_report_trimmed.pdf"
    params:
        pdf="--pdf",
        replace_old="--force" # revisar que no remplaze al anterior
    log:
        "logs/multiqc/multiqc_2.log"
    priority: 6
    shell:
        "multiqc {params.pdf} {params.replace_old} {input} 2>{log}"
