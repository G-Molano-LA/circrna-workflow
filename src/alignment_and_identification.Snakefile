#!/usr/bin/python3

###############################################################################
data= [] # creating an empty list
with open("files.txt") as fp:
    for line in fp:
        data.append(line.strip())
SAMPLES=data
###############################################################################
rule all:
    input:
        "ciri/identification_results",
        "circexplorer2/circularRNA_known.txt"

"""
You need to have at least one rule (target rule) that does not produce any output,
but takes as input all your expected output. Snakemake will take that rule as
first rule, and then checks which rules produce the input of the rule. Then it
checks, which rules produce the inputs for those rules, etc.
"""
## 2. ALIGNMENT AND IDENFITICATION ############################################

rule dw_ref_genome:
    output:
        "raw_data/hg19_ref_genome.fna.gz"
    message:
        " Downloanding reference genome (GRch38)..."
    shell:
        "wget -c https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz "
        "-O {output}"

rule bwa_index:
    input:
        "raw_data/hg19_ref_genome.fna.gz"
    output:
        expand("raw_data/{{genome}}.fna.{ext}", ext=["amb", "ann", "bwt", "pac", "sa"])
    params:
        prefix="{genome}", # no entiendo que he hecho aquí
        algorithm="bwtsw"
    # message:
    #     "Creating a reference genome index.."
    shell:
        "bwa index -a {params.algorithm} -p {params.prefix} {input}"

rule bwa_mem:
    input:
        reads=expand("raw_data/samples/{samples}", samples=SAMPLES), # preguntar esto
        index=expand("raw_data/{genome}.fna.{ext}", genome=["hg19_ref_genome"],
            ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        temp("mapped_reads/{sample}.sam")
    params:
        score="19" # Do not output alignment with score lower than INT.
    # message:
    #     "Executing bwa_mem aligner with {threads} threads on the following files "
    #     "{input}. All alignment with score lower than {params.score} won't be output."
    log:
        "logs/mapped_reads/{sample}_mapped.log"
    shell:
        "bwa mem -T {params.score} {input.index} {input.reads} 1> {output} 2> {log}"

rule ciri2_id:
    input:
        ciri="ciri/CIRI2.pl",
        sam=expand("mapped_reads/{sample}.sam", sample=SAMPLES),
        index=expand("raw_data/{genome}.fna.{ext}", genome=["hg19_ref_genome"],
            ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        "ciri/identification_results"
    # message:
    #     "CIRI2:starting circRNA identification"
    shell:
        "perl {input.ciri} -I {input.sam} -O {output} -F {input.index}"

rule circexplorer2_id:
    input:
        sam=expand("mapped_reads/{sample}.sam", sample=SAMPLES)
    output:
        "circexplorer2/back_spliced_junction.bed"
    log:
        "logs/circexplorer2/parse.log"
    params:
        aligner="BWA"
    # message:
    #     "CircExplorer2: extracting back-spliced exon-exon junction information"
    shell:
        "CIRCexplorer2 parse -t {params.aligner} {input.sam} 2> {log}"

rule circexplorer_annotation:
    input:
        bsj="circexplorer2/back_spliced_junction.bed",
        index=expand("raw_data/{genome}.fna.{ext}", genome=["hg19_ref_genome"],
            ext=["amb", "ann", "bwt", "pac", "sa"]),
        gene="raw_data/GRCh38_latest_genomic.gff.gz"
    output:
        "circexplorer2/circularRNA_known.txt"
    log:
        "logs/circexplorer2/annotate.log"
    # message:
    #     "CircExplorer2: annotating circRNAs with known RefSeq genes"
    shell:
        "CIRCexplorer2 annotate -r {input.gene} -g {input.index} -b {input.bsj} "
        "-o {output} 2> {log}"
