#!/usr/bin/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "ALMOST FINISHED" # requires execution to finish it

# VARIABLES
DATA_TRI  = config["trimming"]["reads"]
TRI_READ1 = expand("{path}/{sample}{ext}", path = config["trimming"]["reads"],
    sample = SAMPLES, ext = config["trimming"]["suffix"][1])
TRI_READ2 = expand("{path}/{sample}{ext}", path = config["trimming"]["reads"],
    sample = SAMPLES, ext = config["trimming"]["suffix"][2])

# TARGET RULE
rule trimming_results:
    input:
        f'{OUTDIR}/quality_control/trimmed/multi_report_trimmed.html'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TRIM_GALORE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rule trimming:
    input:
        read1 = TRI_READ1 if DATA_TRI is not None else RAW_READ1,
        read2 = TRI_READ2 if DATA_TRI is not None else RAW_READ2
    output:
        reads = expand("{outdir}/data/trimmed/{sample}_{replicate}.fq.gz",  outdir = OUTDIR,
            sample = SAMPLES, replicate = ["1_val_1", "2_val_2"]),
        txt   = expand("{outdir}/data/trimmed/{sample}_{replicate}.fastq.gz_trimming_report.txt",
            outdir = OUTDIR, sample = SAMPLES, replicate = [1,2])
    params:
        outdir = f'{OUTDIR}/data/trimmed'
    priority: 8
    shell:
        """
        # Install Trim Galore
        curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz
        tar xvzf trim_galore.tar.gz

        # Running
        ./TrimGalore-0.6.6/trim_galore --paired --output_dir {params.outdir} {input.read1} {input.read2}
        """

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FASTQC_POST-TRIMMING~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule fastqc2:
    input:
        reads = expand("{path}/data/trimmed/{sample}_{replicate}.fq.gz", path = OUTDIR,
            sample = SAMPLES, replicate = ["1_val_1", "2_val_2"])
    output:
        html = expand("{path}/quality_control/trimmed/{sample}_{replicate}_fastqc.html",
            path = OUTDIR, sample = SAMPLES, replicate = ["1_val_1", "2_val_2"]),
        zip  = expand("{path}/quality_control/trimmed/{sample}_{replicate}_fastqc.zip",
            path = OUTDIR, sample = SAMPLES, replicate = ["1_val_1", "2_val_2"])
    params:
        outdir = f'{OUTDIR}/quality_control/trimmed/'
    conda: config["envs"]["quality_control"]
    # message:
    #     "Starting quality analysis control with FASTQC programm on the "
    #     "following files {input.reads}. The number of threads used are {threads}."
    threads: config["quality_control"]["fastqc_threads"]
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
    conda: config["envs"]["quality_control"]
    priority: 6
    shell:
        "multiqc {params.pdf} {params.replace_old} {input} --outdir {params.outdir}"
