#!/usr/bin/python3

## 2. ALIGNMENT AND IDENFITICATION ############################################

rule dw_ref_genome:
    output:
        "raw_data/hg19_ref_genome.fna.gz"
    message: " Downloanding reference genome (GRch38)..."
    shell:
        "wget -c https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz "
        "-O {output}"

rule bwa_index:
    input:
        "raw_data/{genome}.fna.gz"
    output:
        "raw_data/{genome}.amb",
        "raw_data/{genome}.ann",
        "raw_data/{genome}.bwt",
        "raw_data/{genome}.pac",
        "raw_data/{genome}.sa"
    params:
        prefix="{genome}",
        algorithm="bwtsw"
    message: "Creating a reference genome index.."
    shell:
        "bwa index -a {params.algorithm} -p {params.prefix} {input}"

rule bwa_mem:
    input:
        reads=expand("trimmed_data/{sample}_{replicate}.{ext}", replicate=[1,2],
        ext=["fastq.gz", "fastq"]),
        index="raw_data/{genome}"
    output:
        "mapped_reads/{sample}.sam"
    params:
        score="19" # Do not output alignment with score lower than INT.
    message:
        "Executing bwa_mem aligner with {threads} threads on the following files "
        "{input}. All alignment with score lower than {params.score} won't be output."
    log: "mapped_reads/{sample}_mapped.log"
    shell:
        "bwa mem -T {params.score} {input.index} {input.reads} 1> {output} 2> {log}"

rule ciri2_id:
    input:
        ciri="ciri/CIRI2.pl",
        sam="mapped_reads/{sample}.sam",
        index="trimmed_data/hg19_bwa_index.fa"
    output:
        "ciri/identification_results"
    message: "CIRI2:starting circRNA identification"
    shell:
        "perl {input.ciri} -I {input.sam} -O {output} -F {input.index}"

rule circexplorer2_id:
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

rule circexplorer_annotation:
    input:
        bsj="circexplorer2/back_spliced_junction.bed",
        index="trimmed_data/hg19_bwa_index.fa",
        gene="trimmed_data/hg9_gene_all.txt"
    output:
        "circexplorer2/circularRNA_known.txt"
    log:
        "circexplorer2/logs/annotate.log"
    message: "CircExplorer2: annotating circRNAs with known RefSeq genes"
    shell:
        "CIRCexplorer2 annotate -r {input.gene} -g {input.index} -b {input.bsj} "
        "-o {output} 2> {log}"
