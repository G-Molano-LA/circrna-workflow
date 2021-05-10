#!bin/bash/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "IN PROCESS"

DE_DESIGN   = config["DE"]["design"]

rule DE_results:
    input:
        DE_matrix    = f'{OUTDIR}/DE_analysis/circrna_DE.csv',
        volcano_plot = f'{OUTDIR}/DE_analysis/volcano_plot.svg'

rule edgeR:
    input:
        lib_info      = f'{OUTDIR}/DE_analysis/library_info.csv',
        circ_info     = f'{OUTDIR}/DE_analysis/circular_info.csv',
        circ_counts   = f'{OUTDIR}/DE_analysis/circular_count_matrix.csv',
        linear_counts = f'{OUTDIR}/DE_analysis/gene_count_matrix.csv',
        linear_info   = "libs/DE_analysis/linear_info.csv" # comprovar si aquest es genera
    output:
        DE_matrix    = f'{OUTDIR}/DE_analysis/circrna_DE.csv',
        volcano_plot = f'{OUTDIR}/DE_analysis/volcano_plot.svg'
    params:
        design = DE_DESIGN,
        pval   = config["DE"]["pvalue"],
        FC     = config["DE"]["fold_change"],
        outdir = f'{OUTDIR}/DE_analysis',
        script = "src/tools/circ_edge_DE.R",
    conda: config["envs"]["R"]
    shell:
        "Rscript {params.script} -d {params.design}\
                                --lib {input.lib_info}\
                                --circ_info {input.circ_info}\
                                -c {input.circ_counts}\
                                -l {input.linear_counts}\
                                --outdir {params.outdir}"
                                #-g {input.linear_info}"
