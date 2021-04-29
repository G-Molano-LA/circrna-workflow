#!/usr/bin/python3
################################################################################
__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "ALMOST FINISHED" # Falta comprovar funcionamiento con GRCh37
################################################################################
min_version("6.0")

from snakemake.utils import min_version
import subprocess

configfile: "src/utils/config.yaml"

rule results:
    input:
        results = expand("libs/identification/{sample}_coincident_circRNAs.txt",
            sample= config["samples"]),
        diagrams = expand("libs/plots/venn_diagrams/{sample}.png", sample =  config["samples"]),
        ref=expand("data/raw_data/{genome}.fna", genome = config["genome"])

#~~~~~~~~~~~~~~~~~~~~~~~~~ANNOTATION&REFERENCE-FILES~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule get_ref_genome:
    output:
        ref=expand("data/raw_data/{genome}.fna", genome = config["genome"])
    priority: 12
    # message:
    #     "Downloanding"
    run:
        if wildcards.genome == 'GRCh38':
            subprocess.run('wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GO_TO_CURRENT_VERSION/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz \
                            -O data/raw_data/GRCh38.fna.gz \
                            && gunzip data/raw_data/GRCh38.fna.gz \
                            && rm data/raw_data/GRCh38.fna.gz',
                            shell = True, text = True) # UCSC format
        elif wildcards.genome == 'GRCh37':
            subprocess.run('wget -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz \
                            -O data/raw_data/GRCh37.fna.gz \
                            && gunzip data/raw_data/GRCh37.fna.gz \
                            && rm data/raw_data/GRCh37.fna.gz',
                            shell = True, text = True) # ensEMBL format
        else:
            print("Please specify a genome in the config.yaml file. Genome: hg38/GRCh38 or hg37/GRCh37/hg19")


rule get_ref_annotation:
    output:
        expand("data/raw_data/{genome}_ann.gtf", genome = config["genome"])
    priority: 11
    run:
        if wildcards.genome == 'hg38' or 'GRCh38':
            subprocess.run("wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz \
                            -O data/raw_data/GRCh38_ann.gtf.gz \
                            && gunzip data/raw_data/GRCh38_ann.gtf.gz",
                            shell = True, text = True) # UCSC format
        elif wildcards.genome == 'hg37' or 'GRCh37' or 'hg19':
            subprocess.run("wget -c http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz \
                            -O data/raw_data/GRCh37_ann.gtf.gz \
                            && gunzip data/raw_data/GRCh37_ann.gtf.gz",
                            shell = True, text = True) # ensEMBL format


rule refFlat:
    output:
        refFlat=expand("data/raw_data/{genome}_refFlat.txt", genome = config["genome"])
    priority: 10
    run:
        if wildcards.genome == 'hg38' or 'GRCh38':
            subprocess.run("wget -c https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refFlat.txt.gz \
                            -O data/raw_data/GRCh38_refFlat.txt \
                            && gunzip data/raw_data/GRCh38_refFlat.txt",
                            shell = True, text = True)
        elif wildcards.genome == 'hg37' or 'GRCh37' or 'hg19':
            subprocess.run("bash src/utils/gtf_to_refFlat.sh",
                            shell = True, text = True)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~BWA-MEM ALIGNMENT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule bwa_index:
    input:
        ref=expand("data/raw_data/{genome}.fna", genome = config["genome"])
    output:
        index=expand("data/raw_data/bwa/{genome}.fna.{ext}", genome = config["genome"],
            ext=["amb", "ann", "bwt", "pac", "sa"])
    params:
        genome = config["genome"],
        script = "src/utils/bwa_index.sh"
    conda: config["envs"]["bwa"]
    priority: 10
    shell:
        "bash {params.script} {params.genome}"


rule bwa_mem:
    input:
        read1=lambda wildcards: "data/trimmed_data/"+wildcards.sample+config["suffix"][1],
        read2=lambda wildcards: "data/trimmed_data/"+wildcards.sample+config["suffix"][2],
        index=expand("data/raw_data/bwa/{genome}.fna.{ext}", genome = config["genome"],
            ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        sam=temp("data/mapped_data/{sample}.sam")
    params:
        score="19", # Do not output alignment with score lower than INT.
        prefix=expand("data/raw_data/{genome}.fna", genome = config["genome"])
    message:
        "BWA-MEM aligner: Starting the alignment of the reads from {input.read1} & {input.read2}\
         THREADS = {threads}\
         OUTPUT = {output}\
         GENOME = {params.prefix}\
         SCORE = All alignment with score lower than {params.score} won't be output\
         LOG = {log}"
    threads: config["threads"]["bwa"]
    log:
        "logs/mapped_data/{sample}.log"
    conda: config["envs"]["bwa"]
    priority: 7
    shell:
        "bwa mem -T {params.score} -t {threads} {params.prefix} {input.read1} {input.read2} 1> {output.sam} 2> {log}"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CIRI2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule dw_ciri2:
    output:
        "libs/identification/ciri2/CIRI2.pl"
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
        ciri="libs/identification/ciri2/CIRI2.pl",
        sam="data/mapped_data/{sample}.sam",
        ref=expand("data/raw_data/{genome}.fna", genome = config["genome"]) # no acepta archivo comprimido
    output:
        "libs/identification/ciri2/{sample}_results",
    threads: config["threads"]["ciri2"]
    message:
        "CIRI2: Starting circRNA identification in {wildcards.sample}.sam file...\
        REFERENCE FILE = {input.ref}\
        THREADS = {threads} \
        OUTPUT = {output}"
    conda: config["envs"]["ciri2"]
    priority: 5
    shell:
        "perl {input.ciri} -I {input.sam} -O {output} -F {input.ref} -T {threads}"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CIRCEXPLORER2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule circexplorer2_id:
    input:
        sam="data/mapped_data/{sample}.sam"
    output:
        "libs/identification/circexplorer2/{sample}_back_spliced_junction.bed"
    params:
        aligner="BWA"
    message:
        "CIRCEXPLORER2: Extracting back-spliced exon-exon junction information from INPUT = {input} ...\
        ALIGNER = {params.aligner}\
        OUTPUT = {output}"
    conda: config["envs"]["circexplorer2"]
    priority: 6
    shell:
        "CIRCexplorer2 parse -t {params.aligner} --bed={output} {input}"


rule circexplorer2_annotation:
    input:
        bsj="libs/identification/circexplorer2/{sample}_back_spliced_junction.bed",
        ref=expand("data/raw_data/{genome}.fna", genome = config["genome"]),
        refFlat=expand("data/raw_data/{genome}_refFlat.txt", genome = config["genome"])
    output:
        "libs/identification/circexplorer2/{sample}_circularRNA_known.txt"
    message:
        "CIRCEXPLORER2: Annotating circRNAs with known RefSeq genes in {input.bsj}. \
        REFERENCE FILE = {input.ref}\
        ANNOTATING FILE = {input.refFlat} (refFlat format) \
        OUTPUT = {output}"
    conda: config["envs"]["circexplorer2"]
    priority: 5
    shell:
        "CIRCexplorer2 annotate -r {input.refFlat} -g {input.ref} -b {input.bsj} "
        "-o {output}"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OVERLAP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule select_coincidences:
    input:
        ciri2="libs/identification/ciri2/{sample}_results",
        circexplorer2= "libs/identification/circexplorer2/{sample}_circularRNA_known.txt"
    output:
        txt="libs/identification/{sample}_coincident_circRNAs.txt",
        venn="libs/plots/venn_diagrams/{sample}.png"
    params:
        script = "src/utils/select_coincidents.R"
    message:
        "OVERLAP: Selecting circRNAs commonly identified by CIRI2 and CircExplorer2 tools... \
        INPUTS = 1) {input.ciri2} ;  2){input.circexplorer2} \
        OUTPUTS = 1) Venn Diagram: {output.venn} ; 2) Circular matrix: {output.txt} \
        SCRIPT = {params.script}"
    priority: 4
    conda: config["envs"]["R"]
    shell:
        "Rscript {params.script} \
         --ciri2 {input.ciri2} \
         --circexplorer2 {input.circexplorer2}"
