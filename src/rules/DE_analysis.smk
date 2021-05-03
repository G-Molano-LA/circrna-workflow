#!bin/bash/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "IN PROCESS"


rule DE_analysis_results:
    input:

rule edgeR:
    input:
        script="src/main/tools/circ_DE.R",
        lib_info="lib/DE_analysis/library_info.csv",
        circ_info="libs/DE_analysis/circular_info.csv",
        circ_counts="libs/DE_analysis/circular_count_matrix.csv",
        linear_counts="libs/DE_analysis/gene_count_matrix.csv",
        linear_info="libs/DE_analysis/linear_info.csv" # comprovar si aquest es genera
    output:
        "libs/DE_analysis/circrna_DE.csv"
    params:
        design="~group+sex" # això ho ha de poder posar l'investigador
    conda:
        ""
    shell:
        "Rscript {input.script} -d {params.design}\
                                -m {input.lib_info}\
                                --circ_info {input.circ_info}\
                                -c {input.circ_counts}\
                                -l {input.linear_counts}\
                                -g {input.linear_info}"
