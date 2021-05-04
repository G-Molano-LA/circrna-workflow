#!bin/bash/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "IN PROCESS"


rule prep_DE_results:
    input:
        lib_info    = f'{OUTDIR}/DE_analysis/library_info.csv',
        circ_info   = f'{OUTDIR}/DE_analysis/circular_info.csv',
        circ_counts = f'{OUTDIR}/DE_analysis/circular_count_matrix.csv',
        ratio       = f'{OUTDIR}/DE_analysis/junction_ratio.csv',
        gene        = f'{OUTDIR}/DE_analysis/gene_count_matrix.csv',
        transcript  = f'{OUTDIR}/DE_analysis/transcript_count_matrix.csv'


# Preparation of prep_CIRIquant output files & StringTie output files
rule prep_previous_files:
    output:
        circular_file = f'{OUTDIR}/DE_analysis/prep_DE/sample_circ.lst',
        linear_file = f'{OUTDIR}/DE_analysis/prep_DE/sample_gene.lst'
    params:
        dir     = f'{OUTDIR}/DE_analysis/prep_DE',
        outdir  = f'{OUTDIR}/DE_analysis/prep_DE',
        samples = SAMPLES,
        group   = config["prep_DE"]["group"],
        sex     = config["prep_DE"]["sex"],
        script  = "src/main/utils/creating_sample_lst_opt.R"
    priority: 10
    conda: config["envs"]["R"]
    shell:
        "Rscript {params.script} --samples {params.samples}\
            --group {params.group}\
            --sex {params.sex}\
            --dir {params.dir}\
            --outdir {params.outdir} "

rule prep_CIRIquant:
    input:
        circular_file = f'{OUTDIR}/DE_analysis/prep_DE/sample_circ.lst'
    output:
        lib_info    = f'{OUTDIR}/DE_analysis/library_info.csv',
        circ_info   = f'{OUTDIR}/DE_analysis/circular_info.csv',
        circ_counts = f'{OUTDIR}/DE_analysis/circular_count_matrix.csv',
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
