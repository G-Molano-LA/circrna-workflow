#!bin/bash/python3

################################################################################
# Snakefile to explorer the coexpression of circular RNAs with NetMiner tool.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 06-05-2021
# Last modification : 06-05-2021
################################################################################

rule coexpression_results:
    input:
        f'{OUTDIR}/coexpression/final_geneNet.tab'

rule normalization:
    input:
        circ_counts   = f'{OUTDIR}/DE_analysis/circular_count_matrix.csv',
        metadata      = f'{OUTDIR}/DE_analysis/library_info.csv',
        circ_info     = f'{OUTDIR}/DE_analysis/circular_info.csv'
    output:
        expand("{outdir}/data/normalized_counts/{norm}_Count.txt", outdir = OUTDIR,
            norm = ["Raw", "FPKM", "Median", "VST", "TMM", "UQ"])
    params:
        outdir = f'{OUTDIR}/data/normalized_counts/",
        script = 'src/tools/normalization.R'
    conda: config["envs"]["R"]
    shell:
        "Rscript {params.script} --circ_counts {input.circ_counts}\
                --metadata {input.metadata}\
                --circ_info {input.circ_info}\
                --outdir {params.outdir}"
rule netminer:
    input:
        expand("{outdir}/data/normalized_counts/{norm}_Count.txt", outdir = OUTDIR,
            norm = ["Raw", "FPKM", "Median", "VST", "TMM", "UQ"])
    output:
        f'{OUTDIR}/coexpression/final_geneNet.tab'
    params:
        filepath      = config["netminer"]["filepath"],
        percThreshold = config["netminer"]["percThreshold"],
        S1N           = config["netminer"]["S1N"],
        S2N           = config["netminer"]["S2N"],,
        script        = "src/tools/ensemble_method_for_construction_consensus_network.R"
    threads: config["netminer"]["threads"]
    conda: config["envs"]["R"]
    shell:
        "Rscript --vanilla {params.script}  \
        {params.filepath} {threads} {params.percThreshold} {params.S1N} {params.S2N}"
