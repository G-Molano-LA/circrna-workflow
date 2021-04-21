#!bin/bash/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "IN PROCESS"

# Depend on a minimum snakemake version
from snakemake.utils import min_version
min_version("6.0")
###############################################################################
rule results:
    input:

rule previous_files:
    output:
        circular_file = "libs/ciriquant/sample.lst"
        linear_file = "libs/ciriquant/sample_gene.lst"
    priority: 10
    script:
        "src/main/utils/creating_sample_lst.R"

rule prep_CIRIquant:
    input:
        circular_file = "libs/ciriquant/sample.lst"
    output:
        lib_info = "libs/DE_analysis/library_info.csv",
        circ_info = "libs/DE_analysis/circular_info.csv",
        circ_counts = "libs/DE_analysis/circular_count_matrix.csv",
        ratio = "libs/DE_analysis/junction_ratio.csv"
    priority: 9
    shell:
        "prep_CIRIquant -i {input.circular_file} \
                        --lib {output.lib_info} \
                        --circ {output.circ_info} \
                        --bsj {output.circ_counts} \
                        --ratio {output.ratio}"
rule prep_stringtie_output:
    input:
        script="src/main/tools" # ?¿?¿?
        linear_file = "libs/ciriquant/sample_gene.lst"
    output:
        "libs/DE_analysis/gene_count_matrix.csv"
    priority: 9
    shell:
        " {input.script} -i {input.linear_file}"

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
