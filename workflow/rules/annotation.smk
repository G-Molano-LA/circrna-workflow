#!/usr/bin/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "ALMOST FINISHED" # requires execution to finish it

################################################################################
# Snakefile to annotate circRNAs with circBASE ID
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 27-05-2021
# Last modification : 01-06-2021
################################################################################

# VARIABLES
## Config file
LIST   = config["annotation"]["file"]
FIELDS = config["annotation"]["fields"]

# TARGET RULE
rule annotation_results:
    input:
        f'{OUTDIR}/annotation/annotated_circRNA.txt'


rule circbase:
    output:
        f'{OUTDIR}/annotation/hsa_hg19_circRNA.txt'
    params:
        outdir = f'{OUTDIR}/data/trimmed'
    priority: 23
    message: "Downloading circBase database: {output}"
    shell:
        "wget -c http://www.circbase.org/download/hsa_hg19_circRNA.txt -O {output}"

rule chain_file:
    output:
        f'{OUTDIR}/annotation/hg38ToHg19.over.chain'
    params:
        outdir = f'{OUTDIR}/annotation'
    priority: 23
    message:
        "Downloading the chain file provided by UCSC for hg38 to hg19 transformation: {output}"
    shell:
        "wget --timestamping ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz \
            -O {params.outdir}/hg38ToHg19.over.chain.gz\
            && gunzip {params.outdir}/hg38ToHg19.over.chain.gz"

rule annotation:
    input:
        db    = f'{OUTDIR}/annotation/hsa_hg19_circRNA.txt',
        chain = f'{OUTDIR}/annotation/hg38ToHg19.over.chain',
        list  = LIST if LIST is not None else RES_QUANT
    output:
        f'{OUTDIR}/annotation/annotated_circRNA.txt'
    params:
        fields = FIELDS if FIELDS is not None else 'None',
        script = "tools/annotation_db.R"
    conda: config["envs"]["R"]
    priority: 24
    message: "Annotating circRNA ID from circBase database. INPUT: {input.list}; Extra fields selected: {params.fields}"
    shell:
        "Rscript {params.script} --list {input.list} --db {input.db} --chain {input.chain}\
                --out {output} --fields {params.fields}"
