#!bin/python3

################################################################################
# Snakefile to visualize normalized circRNA counts.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 06-05-2021
# Last modification : 07-05-2021
################################################################################
METADATA   = config["visualization"]["metadata"]
COUNT_DATA = config["visualization"]["count_data"]
DESIGN     = config["visualization"]["design"]
FORMAT     = config["visualization"]["out_format"]


rule visualization_results:
    input:
        dendro  = f'{OUTDIR}/visualization/dendrogram.svg',
        boxplot = f'{OUTDIR}/visualization/boxplot.{FORMAT}',
        violin  = f'{OUTDIR}/visualization/violin_plot.{FORMAT}',
        hist    = f'{OUTDIR}/visualization/histogram.{FORMAT}'

rule plots:
    input:
        data = COUNT_DATA if COUNT_DATA != None else CIRC_COUNTS
    output:
        dendro  = f'{OUTDIR}/visualization/dendrogram.svg',
        boxplot = f'{OUTDIR}/visualization/boxplot.{FORMAT}',
        violin  = f'{OUTDIR}/visualization/violin_plot.{FORMAT}',
        hist    = f'{OUTDIR}/visualization/histogram.{FORMAT}'
    params:
        norm      = config["visualization"]["normalized"],
        circ_info = 'None'   if COUNT_DATA != 'None' else CIRC_INFO,
        metadata  = METADATA if COUNT_DATA != 'None' else LIB_INFO,
        design    = DESIGN   if DESIGN     != 'None' else DE_DESIGN,
        output    = FORMAT,
        outdir    = f'{OUTDIR}/visualization',
        script    = "src/tools/visualization.R"
    conda: config["envs"]["R"]
    shell:
        "Rscript {params.script} --data {input.data} --norm {params.norm}\
        --circ_info {params.circ_info} --lib {params.metadata} --output {params.output}\
        --outdir {params.outdir}"
