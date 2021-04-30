#!/usr/bin/python3
################################################################################
__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "ALMOST FINISHED" # Falta comprovar funcionamiento con GRCh37
################################################################################

from snakemake.utils import min_version
import subprocess
min_version("6.0")

# CONFIG FILES
configfile: "src/utils/config.yaml"

# VARIABLES
PATH_genome       = config["path"]["genome_files"]
PATH_trimmed_reads= config["path"]["trimmed_reads"]
OUTDIR            = config["path"]["outdir"]
GENOME            = config["genome"]

# TARGET RULE
rule results:
    input:
        txt  = expand("{path}/identification/overlap/{samples}_common.txt",
            path = OUTDIR, samples = config["samples"]),
        venn = expand("{path}/identification/overlap/{samples}_common.png",
            path = OUTDIR, samples = config["samples"]),
        ref  = f'{PATH_genome}/{GENOME}.fna'

#~~~~~~~~~~~~~~~~~~~~~~~~~ANNOTATION&REFERENCE-FILES~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule get_ref_genome:
    output:
        ref = f'{PATH_genome}/{GENOME}.fna'
    params:
        path = PATH_genome
    priority: 12
    # message:
    #     "Downloanding"
    run:
        if wildcards.genome == 'GRCh38':
            subprocess.run('wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GO_TO_CURRENT_VERSION/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz \
                            -O {params.path}/GRCh38.fna.gz \
                            && gunzip {params.path}/GRCh38.fna.gz \
                            && rm {params.path}/GRCh38.fna.gz',
                            shell = True, text = True) # UCSC format
        elif wildcards.genome == 'GRCh37':
            subprocess.run('wget -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz \
                            -O {params.path}/GRCh37.fna.gz \
                            && gunzip {params.path}/GRCh37.fna.gz \
                            && rm {params.path}/GRCh37.fna.gz',
                            shell = True, text = True) # ensEMBL format
        else:
            print("Please specify a genome in the config.yaml file. Genome: hg38/GRCh38 or hg37/GRCh37/hg19")


rule get_ref_annotation:
    output:
        f'{PATH_genome}/{GENOME}_ann.gtf'
    params:
        path = PATH_genome
    priority: 11
    run:
        if wildcards.genome == 'hg38' or 'GRCh38':
            subprocess.run("wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz \
                            -O {params.path}/GRCh38_ann.gtf.gz \
                            && gunzip {params.path}/GRCh38_ann.gtf.gz",
                            shell = True, text = True) # UCSC format
        elif wildcards.genome == 'hg37' or 'GRCh37' or 'hg19':
            subprocess.run("wget -c http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz \
                            -O {params.path}/GRCh37_ann.gtf.gz \
                            && gunzip {params.path}/GRCh37_ann.gtf.gz",
                            shell = True, text = True) # ensEMBL format


rule refFlat:
    output:
        refFlat = f'{PATH_genome}/{GENOME}.refFlat.txt'
    params:
        path = PATH_genome
    priority: 10
    run:
        if wildcards.genome == 'hg38' or 'GRCh38':
            subprocess.run("wget -c https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refFlat.txt.gz \
                            -O {params.path}/GRCh38_refFlat.txt \
                            && gunzip {params.path}/GRCh38_refFlat.txt",
                            shell = True, text = True)
        elif wildcards.genome == 'hg37' or 'GRCh37' or 'hg19':
            subprocess.run("bash src/utils/gtf_to_refFlat.sh {params.path}",
                            shell = True, text = True)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~BWA-MEM ALIGNMENT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule bwa_index:
    input:
        ref = f'{PATH_genome}/{GENOME}.fna'
    output:
        index=expand("{path}/bwa/{genome}.fna.{ext}", path = PATH_genome,
            genome = GENOME, ext=["amb", "ann", "bwt", "pac", "sa"])
    params:
        genome = GENOME,
        script = "src/utils/bwa_index.sh",
        path = PATH_genome
    conda: config["envs"]["bwa"]
    priority: 10
    shell:
        "bash {params.script} {params.genome} {params.path}"


rule bwa_mem:
    input:
        read1 = lambda wildcards: f'{PATH_trimmed_reads}/{wildcards.sample}{config["suffix"]["trimmed"][1]}',
        read2 = lambda wildcards: f'{PATH_trimmed_reads}/{wildcards.sample}{config["suffix"]["trimmed"][2]}',
        index = expand("{path}/bwa/{genome}.fna.{ext}", genome = GENOME,
            path = PATH_genome, ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        sam    = temp(f'{OUTDIR}/data/mapped_data/{{sample}}.sam')
    params:
        score  = "19", # Do not output alignment with score lower than INT.
        prefix = f'{PATH_genome}/{GENOME}.fna'
    message:
        "BWA-MEM aligner: Starting the alignment of the reads from {input.read1} & {input.read2}\
         THREADS = {threads}\
         OUTPUT  = {output}\
         GENOME  = {params.prefix}\
         SCORE   = All alignment with score lower than {params.score} won't be output\
         LOG     = {log}"
    threads: config["threads"]["bwa"]
    log:
        f'{OUTDIR}logs/mapped_data/{{sample}}.log'
    conda: config["envs"]["bwa"]
    priority: 7
    shell:
        "bwa mem -T {params.score} -t {threads} {params.prefix} {input.read1} {input.read2} 1> {output.sam} 2> {log}"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CIRI2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule dw_ciri2:
    output:
        "src/tools/CIRI2.pl"
    message:
        "Downloanding alignment tool CIRI2 from official website...\
        SCRIPT = {params.script}"
    params:
        script = "src/utils/install_ciri2.sh"
    priority: 6
    shell:
        "bash {params.script}"


rule ciri2:
    input:
        ciri = "src/tools/CIRI2.pl",
        sam  = f'{OUTDIR}/data/mapped_data/{{sample}}.sam',
        ref  =  f'{PATH_genome}/{GENOME}.fna' # no acepta archivo comprimido
    output:
        f'{OUTDIR}/identification/{{sample}}_ciri2.txt'
    threads: config["threads"]["ciri2"]
    message:
        "CIRI2: Starting circRNA identification in {wildcards.sample}.sam file...\
        REFERENCE FILE = {input.ref}\
        THREADS        = {threads} \
        OUTPUT         = {output}"
    conda: config["envs"]["ciri2"]
    priority: 5
    shell:
        "perl {input.ciri} -I {input.sam} -O {output} -F {input.ref} -T {threads}"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CIRCEXPLORER2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule circexplorer2_id:
    input:
        sam = f'{OUTDIR}/data/mapped_data/{{sample}}.sam',
    output:
        temp(f'{OUTDIR}/identification/{{sample}}_back_spliced_junction.bed')
    params:
        aligner="BWA"
    message:
        "CIRCEXPLORER2: Extracting back-spliced exon-exon junction information from INPUT = {input} ...\
        ALIGNER = {params.aligner}\
        OUTPUT  = {output}"
    conda: config["envs"]["circexplorer2"]
    priority: 6
    shell:
        "CIRCexplorer2 parse -t {params.aligner} --bed={output} {input}"


rule circexplorer2_annotation:
    input:
        bsj     = f'{OUTDIR}/identification/{{sample}}_back_spliced_junction.bed',
        ref     = f'{PATH_genome}/{GENOME}.fna',
        refFlat = f'{PATH_genome}/{GENOME}.refFlat.txt'
    output:
        f'{OUTDIR}/identification/{{sample}}_circexp2.txt'
    message:
        "CIRCEXPLORER2: Annotating circRNAs with known RefSeq genes in {input.bsj}. \
        REFERENCE FILE  = {input.ref}\
        ANNOTATING FILE = {input.refFlat} (refFlat format) \
        OUTPUT          = {output}"
    conda: config["envs"]["circexplorer2"]
    priority: 5
    shell:
        "CIRCexplorer2 annotate -r {input.refFlat} -g {input.ref} -b {input.bsj} "
        "-o {output}"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OVERLAP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule select_coincidences:
    input:
        ciri2         = f'{OUTDIR}/identification/{{sample}}_ciri2.txt',
        circexplorer2 = f'{OUTDIR}/identification/{{sample}}_circexp2.txt'
    output:
        txt  = f'{OUTDIR}/identification/overlap/{{sample}}_common.txt',
        venn = f'{OUTDIR}/identification/overlap/{{sample}}_common.png'
    params:
        script = "src/utils/select_coincidents.R",
    message:
        "OVERLAP: Selecting circRNAs commonly identified by CIRI2 and CircExplorer2 tools... \
        INPUTS  = 1) {input.ciri2} ;  2){input.circexplorer2} \
        OUTPUTS = 1) Venn Diagram: {output.venn} ; 2) Circular matrix: {output.txt} \
        SCRIPT  = {params.script}"
    priority: 4
    conda: config["envs"]["R"]
    shell:
        "Rscript {params.script} \
         --ciri2 {input.ciri2} \
         --circexplorer2 {input.circexplorer2}"
