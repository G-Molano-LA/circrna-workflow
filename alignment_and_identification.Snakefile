#!/usr/bin/python3

## 2. ALIGNMENT AND IDENFITICATION ############################################

rule bwa_index
    input:
        reads=expand("raw_data/{sample}_{replicate}.{ext}", replicate=[1,2],
        ext=["fastq.gz", "fastq"])
    output:
        "raw_data/hg19_bwa_index.fa"
    message: "Creating a reference genome index"
    shell:
        "bwa index -a bwtsw hg19_bwa_index.fa"

rule bwa_mem:
    input:
        reads=expand("raw_data/{sample}_{replicate}.{ext}", replicate=[1,2],
        ext=["fastq.gz", "fastq"])
        index="raw_data/hg19_bwa_index.fa"
    output:
        "mapped_reads/{sample}.sam"
    params:
        score="19" # Do not output alignment with score lower than INT.
    threads: 8
    message:
        "Executing bwa_mem aligner with {threads} threads on the following files "
        "{input}. All alignment with score lower than {params.score} won't be output."
    log: "mapped_reads/{sample}_mapped.log"
    shell:
        "bwa mem -T {params.score} -t {threads} {input.index} {input.reads} 1> {output} 2> {log}"

rule ciri2_id
    input:
        ciri="ciri/CIRI2.pl"
        sam="mapped_reads/{sample}.sam"
        index="raw_data/hg19_bwa_index.fa"
    output:
        "ciri/identification_results"
    message: "CIRI2:starting circRNA identification"
    shell:
        "perl {input.ciri} -I {input.sam} -O {output} -F {input.index}"

rule circexplorer2_id
    input:
        sam="mapped_reads/{sample}.sam"
    output:
        "circexplorer2/back_spliced_junction.bed"
    log:
        "circexplorer2/logs/parse.log"
    params:
        aligner="BWA"
    message: "CircExplorer2: extracting back-spliced exon-exon junction information"
    shell:
        "CIRCexplorer2 parse -t {params.aligner} {input.sam} 2> {log}"

rule circexplorer_annotation
    input:
        bsj="circexplorer2/back_spliced_junction.bed"
        index="raw_data/hg19_bwa_index.fa"
        gene="raw_data/hg9_gene_all.txt"
    output:
        "circexplorer2/circularRNA_known.txt"
    log:
        "circexplorer2/logs/annotate.log"
    message: "CircExplorer2: annotating circRNAs with known RefSeq genes"
    shell:
        "CIRCexplorer2 annotate -r {input.gene} -g {input.index} -b {input.bsj} "
        "-o {output} 2> {log}"
