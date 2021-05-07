#!bin/bash/python3

################################################################################
# Snakefile to visualize normalized circRNA counts.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 06-05-2021
# Last modification : 06-05-2021
################################################################################
METADATA   = config["visualization"]["metadata"]
COUNT_data = config["visualization"]["count_data"]

rule visualization_results:
    input:
        directory(f'{OUTDIR}/visualization')

rule plots:
    input:
        data = COUNT_data if COUNT_data != None else f'{OUTDIR}/DE_analysis/circular_count_matrix.csv'
    output:
        directory(f'{OUTDIR}/visualization')
    params:
        norm      = config["visualization"]["normalized"],
        circ_info = 'None' if COUNT_data != 'None' else f'{OUTDIR}/DE_analysis/circular_info.csv',
        metadata  = METADATA if COUNT_data != 'None' else f'{OUTDIR}/DE_analysis/library_info.csv',
        output    = config["visualization"]["out_format"],
        script    = "src/tools/visualization.R"
    conda: config["envs"]["R"]
    shell:
        "Rscript {params.script} --data {input.data} --norm {params.norm}\
        --circ_info {params.circ_info} --lib {params.metadata} --output {params.output}\
        --outdir {output}"
