#!bin/python3

################################################################################
# Snakefile to visualize normalized circRNA counts.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 06-05-2021
# Last modification : 01-06-2021
################################################################################
METADATA_VI = config["visualization"]["metadata"]
DATA_VI     = config["visualization"]["count_data"]
DESIGN_VI   = config["visualization"]["design"]
FORMAT_VI   = config["visualization"]["out_format"]


rule visualization_results:
    input:
        dendro  = f'{OUTDIR}/visualization/dendrogram.svg',
        boxplot = f'{OUTDIR}/visualization/boxplot.{FORMAT_VI}',
        violin  = f'{OUTDIR}/visualization/violin_plot.{FORMAT_VI}',
        hist    = f'{OUTDIR}/visualization/histogram.{FORMAT_VI}'

rule plots:
    input:
        data = CIRC_COUNTS if DATA_VI is None else DATA_VI
    output:
        dendro  = f'{OUTDIR}/visualization/dendrogram.svg',
        boxplot = f'{OUTDIR}/visualization/boxplot.{FORMAT_VI}',
        violin  = f'{OUTDIR}/visualization/violin_plot.{FORMAT_VI}',
        hist    = f'{OUTDIR}/visualization/histogram.{FORMAT_VI}'
    params:
        norm      = config["visualization"]["normalized"],
        metadata  = LIB_INFO if METADATA_VI is None else METADATA_VI,
        separator = config["visualization"]["metadata_sep"],
        design    = DE_DESIGN if DESIGN_VI is None else DESIGN_VI,
        output    = FORMAT_VI,
        units     = config["visualization"]["units"],
        width     = config["visualization"]["width"],
        height    = config["visualization"]["height"],
        outdir    = f'{OUTDIR}/visualization',
        script    = "tools/visualization.R"
    conda: config["envs"]["R"]
    priority: 29
    shell:
        "Rscript {params.script} --data {input.data} --norm {params.norm}\
        --metadata {params.metadata} --sep {params.separator}\
        --output {params.output}\
        --units {params.units} --height {params.height} --width {params.width}\
        --outdir {params.outdir}"
