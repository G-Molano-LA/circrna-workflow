#!bin/bash/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "IN PROCESS"


rule prep_results:
    input:

rule previous_files:
    output:
        circular_file = f'{OUTDIR}/DE_analysis/prep_DE/sample_circ.lst',
        linear_file = f'{OUTDIR}/DE_analysis/prep_DE/sample_gene.lst'
    params:
        dir     = f'{OUTDIR}/DE_analysis/prep_DE',
        outdir  = f'{OUTDIR}/DE_analysis/prep_DE',
        samples = SAMPLES
        group   = config["prep_DE"]["group"]
        sex     = config["prep_DE"]["sex"]
    priority: 10
    script:
        "src/main/utils/creating_sample_lst_opt.R  --samples {params.samples}\
            --dir {params.dir}\
            --outdir {params.outdir}\
            --group {params.group}\
            --sex {params.sex}"

rule prep_CIRIquant:
    input:
        circular_file = f'{OUTDIR}/ciriquant/sample.lst'
    output:
        lib_info    = f'{OUTDIR}/DE_analysis/library_info.csv',
        circ_info   = f'{OUTDIR}/DE_analysis/circular_info.csv',
        circ_counts = f'{OUTDIR}/DE_analysis/circular_count_matrix.csv',
        ratio       = f'{OUTDIR}/DE_analysis/junction_ratio.csv'
    priority: 9
    shell:
        "prep_CIRIquant -i {input.circular_file} \
                        --lib {output.lib_info} \
                        --circ {output.circ_info} \
                        --bsj {output.circ_counts} \
                        --ratio {output.ratio}"
rule prep_stringtie_output:
    input:
        script="src/main/tools", # ?¿?¿?
        linear_file = "libs/ciriquant/sample_gene.lst"
    output:
        "libs/DE_analysis/gene_count_matrix.csv"
    priority: 9
    shell:
        " {input.script} -i {input.linear_file}"
