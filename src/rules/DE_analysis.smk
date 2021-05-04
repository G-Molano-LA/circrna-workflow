#!bin/bash/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "IN PROCESS"


rule DE_results:
    input:
        f'{OUTDIR}/DE_analysis/circrna_DE.csv'

rule edgeR:
    input:
        lib_info      = f'{OUTDIR}/DE_analysis/library_info.csv',
        circ_info     = f'{OUTDIR}/DE_analysis/circular_info.csv',
        circ_counts   = f'{OUTDIR}/DE_analysis/circular_count_matrix.csv',
        linear_counts = f'{OUTDIR}/DE_analysis/gene_count_matrix.csv',
        linear_info   = "libs/DE_analysis/linear_info.csv" # comprovar si aquest es genera
    output:
        f'{OUTDIR}/DE_analysis/circrna_DE.csv'
    params:
        script = "src/tools/circ_edge_DE.R",
        design = config["DE"]["design"],
    conda: config["envs"]["R"]
    shell:
        "Rscript {params.script} -d {params.design}\
                                -m {input.lib_info}\
                                --circ_info {input.circ_info}\
                                -c {input.circ_counts}\
                                -l {input.linear_counts}\
                                -g {input.linear_info}"
