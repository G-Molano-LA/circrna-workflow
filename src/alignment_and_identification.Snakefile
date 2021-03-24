#!/usr/bin/python3

from snakemake.utils import min_version
min_version("6.0")

__author__ = "G. Molano, LA (gonmola@hotmail.es)"

# 1. DATA######################################################################
data= [] # creating an empty list
with open("data/seqs.txt") as fp:
    for line in fp:
        data.append(line.rstrip('_[1,2]_val_[1,2].fq.gz\n'))
data=list(dict.fromkeys(data)) # removing duplicates
SAMPLES=data
###############################################################################

rule results:
    input:
        "libs/ciri2/ciri2_merged_results",
        "libs/circexplorer2/circexplorer2_merged_results.txt"



"""
You need to have at least one rule (target rule) that does not produce any output,
but takes as input all your expected output. Snakemake will take that rule as
first rule, and then checks which rules produce the input of the rule. Then it
checks, which rules produce the inputs for those rules, etc.
"""
# 2. ALIGNMENT AND IDENFITICATION #############################################

rule dw_ref_genome:
    output:
        "data/raw_data/hg19_ref.fna"
    message:
        " Downloanding reference genome (GRch38)..."
    shell:
        "wget -c https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz "
        "-O data/raw_data/hg19_ref.fna.gz && gunzip data/raw_data/hg19_ref.fna.gz"

rule dw_ref_annotation:
    output:
        "data/raw_data/hg19_ann.gff"
    message:
        "Downloanding annotation reference ..."
    shell:
        "wget -c https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz "
        "-O data/raw_data/hg19_ann.gff.gz && gunzip data/raw_data/hg19_ann.gff.gz "

rule bwa_index:
    input:
        "data/raw_data/hg19_ref.fna"
    output:
        expand("data/raw_data/{genome}.fna.{ext}", genome=["hg19_ref"],
            ext=["amb", "ann", "bwt", "pac", "sa"])
    params:
        prefix="hg19_ref",
        algorithm="bwtsw"
    message:
        "Creating a reference genome index from {params.prefix} file"
    log:
         "logs/bwa_index/hg19_ref.log"
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa index -a {params.algorithm} -p {params.prefix} {input} 1> {output} 2> {log}"

rule bwa_mem:
    input:
        reads1=lambda wildcards: "data/trimmed_data/"+wildcards.sample+"_1_val_1.fq.gz",
        reads2=lambda wildcards: "data/trimmed_data/"+wildcards.sample+"_2_val_2.fq.gz"
    output:
        sam="data/mapped_data/{sample}.sam"
    params:
        score="19", # Do not output alignment with score lower than INT.
        prefix="data/raw_data/hg19_ref"
    message:
        "Executing bwa_mem aligner with {threads} threads on the following files "
        "{input}. All alignment with score lower than {params.score} won't be output."
    threads: 10
    log:
        "logs/mapped_data/{sample}.log"
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa mem -T {params.score} -t {threads} {params.prefix} {input.reads1} {input.reads2} 1> {output.sam} 2> {log}"

rule dw_ciri2:
    output:
        "libs/ciri2/CIRI2.pl"
    message:
        " Downloanding alignment tool CIRI2..."
    shell:
        ". src/install_ciri2.sh"

rule ciri2_id:
    input:
        ciri="libs/ciri2/CIRI2.pl",
        sam="data/mapped_data/{sample}.sam",
        ref="data/raw_data/hg19_ref.fna" # no acepta archivo comprimido
    output:
        "libs/ciri2/{sample}_results"
    message:
        "CIRI2:starting circRNA identification in {sample} file"
    log:
        "logs/ciri2/{sample}.log"
    conda:
        "envs/ciri2.yaml"
    shell:
        "perl {input.ciri} -I {input.sam} -O {output} -F {input.ref} 2> {log}"

rule ciri2_merge:
    input:
        expand("libs/ciri2/{sample}_results", sample=SAMPLES)
    output:
        "libs/ciri2/ciri2_merged_results"
    message:
        "Merging results from ciri2 identification ({input}) in one file: {output}"
    params:
        algorithm="ciri2"
    shell:
        "python3 src/circ_merge.py -f {input} -a {params.algorithm} > {output}"


rule circexplorer2_id:
    input:
        sam="data/mapped_data/{sample}.sam"
    output:
        "libs/circexplorer2/{sample}_back_spliced_junction.bed"
    log:
        "logs/circexplorer2/{sample}_parse.log"
    params:
        aligner="BWA"
    message:
        "CircExplorer2: extracting back-spliced exon-exon junction information from {input.sam}"
    conda:
        "envs/circexplorer2.yaml"
    shell:
        "CIRCexplorer2 parse -t {params.aligner} {input} --bed={output} 2> {log}"

rule circexplorer2_annotation:
    input:
        bsj="libs/circexplorer2/{sample}_back_spliced_junction.bed",
        ref="data/raw_data/hg19_ref.fna",
        gene="data/raw_data/hg19_ann.gff"
    output:
        "libs/circexplorer2/{sample}_circularRNA_known.txt"
    log:
        "logs/circexplorer2/{sample}_annotate.log"
    message:
        "CircExplorer2: annotating circRNAs with known RefSeq genes in {input.bsj}. ",
        "The output file {output} will be generated."
    conda:
        "envs/circexplorer2.yaml"
    shell:
        "CIRCexplorer2 annotate -r {input.gene} -g {input.ref} -b {input.bsj} "
        "-o {output} 2> {log}"

rule circexplorer2_merge:
    input:
        expand("libs/circexplorer2/{sample}_circularRNA_known.txt", sample=SAMPLES)
    output:
        "libs/circexplorer2/circexplorer2_merged_results.txt"
    message:
        "Merging results from circexplorer2 identification ({input}) in one file: {output}"
    params:
        algorithm="circexplorer2"
    shell:
        "python3 src/circ_merge.py -f {input} -a {params.algorithm} > {output}"
