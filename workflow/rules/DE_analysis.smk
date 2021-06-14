#!bin/bash/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__  = "IN PROCESS"
################################################################################
# Snakefile to visualize to prepare output CIRIquant files to DE analysis.
# Generation of all sample count matrixs for linear and circular RNAs, respectively.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 03-05-2021
# Last modification : 01-06-2021
################################################################################

# VARIABLES
## Generated
LIB_INFO    = f'{OUTDIR}/DE_analysis/library_info.csv'
CIRC_INFO   = f'{OUTDIR}/DE_analysis/circular_info.csv'
CIRC_COUNTS = f'{OUTDIR}/DE_analysis/circular_count_matrix.csv'
## Config file
METADATA_DE    = config["DE"]["metadata"]
DE_DESIGN      = config["DE"]["design"]


# TARGET RULE
rule DE_results:
    input:
        DE_matrix    = f'{OUTDIR}/DE_analysis/circrna_DE.csv',
        volcano_plot = f'{OUTDIR}/DE_analysis/volcano_plot.svg'

#~~~~~~~~~~~~~~~~~~~~~~~~~~PREPARATION_FILES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Preparation of prep_CIRIquant output files & StringTie output files
rule prep_previous_files:
    output:
        circular_file = f'{OUTDIR}/DE_analysis/prep_DE/sample_circ.lst',
        linear_file   = f'{OUTDIR}/DE_analysis/prep_DE/sample_gene.lst'
    params:
        dir      = f'{OUTDIR}/ciriquant', # directory of ciriquant results
        outdir   = f'{OUTDIR}/DE_analysis/prep_DE',
        metadata = METADATA_DE,
        separator= config["DE"]["metadata_sep"],
        script   = "utils/creating_sample_lst_opt.R"
    priority: 25
    conda: config["envs"]["R"]
    message: "Preparation of "
    shell:
        "Rscript {params.script} --metadata {params.metadata}\
            --sep {params.separator}\
            --dir {params.dir}\
            --outdir {params.outdir}"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~CIRCULAR_RNAs_COUNTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule prep_CIRIquant:
    input:
        circular_file = f'{OUTDIR}/DE_analysis/prep_DE/sample_circ.lst'
    output:
        lib_info    = LIB_INFO,
        circ_info   = CIRC_INFO,
        circ_counts = CIRC_COUNTS,
        ratio       = f'{OUTDIR}/DE_analysis/junction_ratio.csv'
    priority: 26
    conda: config["envs"]["ciriquant"]
    shell:
        "prep_CIRIquant -i {input.circular_file} \
                        --lib {output.lib_info} \
                        --circ {output.circ_info} \
                        --bsj {output.circ_counts} \
                        --ratio {output.ratio}"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~LINEAR_RNAs_COUNTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule StringTie:
    input:
        linear_file = f'{OUTDIR}/DE_analysis/prep_DE/sample_gene.lst',
        script      = "utils/prepDE.py"
    output:
        gene       = f'{OUTDIR}/DE_analysis/gene_count_matrix.csv',
        transcript = f'{OUTDIR}/DE_analysis/transcript_count_matrix.csv'
    priority: 27
    conda: config["envs"]["ciriquant"]
    shell:
        "{input.script} -i {input.linear_file} -g {output.gene} -t {output.transcript}"

#~~~~~~~~~~~~~~~~~~~~~~~~~DIFFERENTIAL_EXPRESSION_ANALYSIS~~~~~~~~~~~~~~~~~~~~~~
rule edgeR:
    input:
        lib_info      = LIB_INFO,
        circ_info     = CIRC_INFO,
        circ_counts   = CIRC_COUNTS,
        linear_counts = f'{OUTDIR}/DE_analysis/gene_count_matrix.csv'
    output:
        DE_matrix    = f'{OUTDIR}/DE_analysis/circrna_DE.csv',
        volcano_plot = f'{OUTDIR}/DE_analysis/volcano_plot.svg',
        heatmap      = f'{OUTDIR}/DE_analysis/heatmap.svg'
    params:
        design = DE_DESIGN,
        pval   = config["DE"]["pvalue"],
        FC     = config["DE"]["fold_change"],
        outdir = f'{OUTDIR}/DE_analysis',
        script = "tools/circ_edge_DE.R",
    conda: config["envs"]["R"]
    priority: 28
    shell:
        "Rscript {params.script} --design {params.design}\
                                --metadata {input.lib_info}\
                                --circ_info {input.circ_info}\
                                --circ_counts {input.circ_counts}\
                                --linear_counts {input.linear_counts}\
                                --pval {params.pval}\
                                --fc {params.FC}\
                                --outdir {params.outdir}"
