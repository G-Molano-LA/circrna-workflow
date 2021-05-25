#!bin/bash/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__  = "IN PROCESS"

METADATA_DE    = config["DE"]["metadata"]

LIB_INFO    = f'{OUTDIR}/DE_analysis/library_info.csv'
CIRC_INFO   = f'{OUTDIR}/DE_analysis/circular_info.csv'
CIRC_COUNTS = f'{OUTDIR}/DE_analysis/circular_count_matrix.csv'

rule prep_DE_results:
    input:
        lib_info    = LIB_INFO,
        circ_info   = CIRC_INFO,
        circ_counts = CIRC_COUNTS,
        ratio       = f'{OUTDIR}/DE_analysis/junction_ratio.csv',
        gene        = f'{OUTDIR}/DE_analysis/gene_count_matrix.csv',
        transcript  = f'{OUTDIR}/DE_analysis/transcript_count_matrix.csv'


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
        script   = "src/utils/creating_sample_lst_opt.R"
    priority: 10
    conda: config["envs"]["R"]
    shell:
        "Rscript {params.script} --metadata {params.metadata}\
            --sep {params.separator}\
            --dir {params.dir}\
            --outdir {params.outdir}"

rule prep_CIRIquant:
    input:
        circular_file = f'{OUTDIR}/DE_analysis/prep_DE/sample_circ.lst'
    output:
        lib_info    = LIB_INFO,
        circ_info   = CIRC_INFO,
        circ_counts = CIRC_COUNTS,
        ratio       = f'{OUTDIR}/DE_analysis/junction_ratio.csv'
    priority: 9
    conda: config["envs"]["ciriquant"]
    shell:
        "prep_CIRIquant -i {input.circular_file} \
                        --lib {output.lib_info} \
                        --circ {output.circ_info} \
                        --bsj {output.circ_counts} \
                        --ratio {output.ratio}"

rule StringTie:
    input:
        linear_file = f'{OUTDIR}/DE_analysis/prep_DE/sample_gene.lst',
        script      = "src/utils/prepDE.py"
    output:
        gene       = f'{OUTDIR}/DE_analysis/gene_count_matrix.csv',
        transcript = f'{OUTDIR}/DE_analysis/transcript_count_matrix.csv'
    priority: 9
    conda: config["envs"]["ciriquant"]
    shell:
        "{input.script} -i {input.linear_file} -g {output.gene} -t {output.transcript}"
