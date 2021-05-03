#!bin/bash/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "NOT STARTED"

# Depend on a minimum snakemake version
from snakemake.utils import min_version
min_version("6.0")
###############################################################################
rule visualization_results:
    input:


rule plots:
    input:
        script="src/circ_DE.R",
        data="lib/DE_analysis/circrna_DE.csv",
    output:
        "libs/DE_analysis/circrna_DE.{params.output}"
    params:
        output="pdf",
        pval="",
        FC=""
    conda:
        ""
    shell:
        "Rscript {input.script} \
        -d {input.data}\
        -p {params.pval}\
        -FC {params.FC}\
        -o {params.output}"
