#!bin/bash/python3

################################################################################
# Snakefile to explorer the coexpression of circular RNAs with NetMiner tool.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 06-05-2021
# Last modification : 01-06-2021
################################################################################
DATA_CO     = config["netminer"]["count_data"]
METADATA_CO = config["netminer"]["metadata"]
DESIGN_CO   = config["netminer"]["design"]
GENE_LENGTH  = config["netminer"]["gene_length"]

rule coexpression_results:
    input:
        f'{OUTDIR}/coexpression/final_geneNet.tab'

rule normalization:
    input:
        circ_counts   = CIRC_COUNTS if DATA_CO is None else DATA_CO ,
        metadata      = METADATA_DE if METADATA_CO is None else METADATA_CO
    output:
        expand("{outdir}/data/normalized_counts/{norm}_Count.txt", outdir = OUTDIR,
            norm = ["Raw", "FPKM", "Median", "VST", "TMM", "UQ"])
    params:
        design    = DE_DESIGN if DESIGN_CO is None else DESIGN_CO,
        circ_info = CIRC_INFO if DATA_CO is None else GENE_LENGTH,
        separator = config["netminer"]["metadata_sep"],
        outdir    = f'{OUTDIR}/data/normalized_counts/',
        script    = 'src/tools/normalization.R'
    conda: config["envs"]["R"]
    priority: 73
    shell:
        "Rscript {params.script} --circ_counts {input.circ_counts}\
                --metadata {input.metadata}\
                --sep {params.separator}\
                --design {params.design}\
                --circ_info {params.circ_info}\
                --outdir {params.outdir}"
rule netminer:
    input:
        expand("{outdir}/data/normalized_counts/{norm}_Count.txt", outdir = OUTDIR,
            norm = ["Raw", "FPKM", "Median", "VST", "TMM", "UQ"])
    output:
        f'{OUTDIR}/coexpression/final_geneNet.tab'
    params:
        filepath      = f'{OUTDIR}/data/normalized_counts',
        percThreshold = config["netminer"]["percThreshold"],
        S1N           = config["netminer"]["S1N"],
        S2N           = config["netminer"]["S2N"],
        script        = "src/tools/ensemble_method_for_construction_consensus_network.R"
    threads: config["netminer"]["threads"]
    conda: config["envs"]["R"]
    priority: 72
    shell:
        "Rscript --vanilla {params.script}  \
        {params.filepath} {threads} {params.percThreshold} {params.S1N} {params.S2N}"
