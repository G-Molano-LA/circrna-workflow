#!bin/bash/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "IN PROCESS"

# Depend on a minimum snakemake version
from snakemake.utils import min_version
min_version("6.0")
###############################################################################
rule results:
    input:


rule edgeR:
    input:
        script="src/tools/circ_DE.R",
        metadata="lib/DE_analysis/metadata_samples.csv",
        circ_counts="libs/DE_analysis/circular_counts_sub.csv",
        linear_counts="libs/DE_analysis/linear_counts_sub.csv",
        linear_info="libs/DE_analysis/linear_info.csv",
        circ_info="libs/DE_analysis/circular_info.csv"
    output:
        "libs/DE_analysis/circrna_DE.csv"
    params:
        design="~group+sex"
    conda:
        ""
    shell:
        "Rscript {input.script} \
        -d {params.design}\
        -m {input.metadata}\
        -c {input.circ_counts}\
        -l {input.linear_counts}\
        -g {input.linear_info}\
        --circ_info {input.circ_info}"
