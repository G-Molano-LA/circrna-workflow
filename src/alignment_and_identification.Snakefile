#!/usr/bin/python3


###############################################################################
data= [] # creating an empty list
with open("files.txt") as fp:
    for line in fp:
        data.append(line.rstrip('_[1,2].fastq.gz\n'))
data=list(dict.fromkeys(data)) # removing duplicates
SAMPLES=data
###############################################################################
rule results:
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
    log:
         "logs/bwa_index/{genome}.log"
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa index -a {params.algorithm} -p {params.prefix} {input} 2> {log}"

rule bwa_mem:
    input:
        reads=expand("trimmed_data/{sample}_{replicate}.fq.gz",
             sample=SAMPLES, replicate=["1_val_1", "2_val_2"]),
        index=expand("raw_data/{genome}.fna.{ext}", genome=["hg19_ref_genome"],
            ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        temp(expand("mapped_reads/{sample}_{replicate}.sam",
            sample=SAMPLES, replicate=[1,2])
            )
    params:
        score="19" # Do not output alignment with score lower than INT.
    # message:
    #     "Executing bwa_mem aligner with {threads} threads on the following files "
    #     "{input}. All alignment with score lower than {params.score} won't be output."
    log:
        expand("logs/mapped_reads/{sample}_{replicate}_mapped.log",
            sample=SAMPLES, replicate=[1,2])
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa mem -T {params.score} {input.index} {input.reads} 1> {output} 2> {log}"

rule dw_ciri2:
    output:
        "ciri2/CIRI2.pl"
    # message:
    #     " Downloanding alignment tool CIRI2..."
    shell:
        ". src/install_ciri2.sh"

rule ciri2_id:
    input:
        ciri="ciri2/CIRI2.pl",
        sam=expand("mapped_reads/{sample}_{replicate}.sam",
            sample=SAMPLES, replicate=[1,2]),
        index=expand("raw_data/{genome}.fna.{ext}", genome=["hg19_ref_genome"],
            ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        "ciri/identification_results"
    # message:
    #     "CIRI2:starting circRNA identification"
    log:
        "logs/ciri2/identification.log"
    conda:
        "envs/ciri2.yaml"
    shell:
        "perl {input.ciri} -I {input.sam} -O {output} -F {input.index} 2> {log}"

rule circexplorer2_id:
    input:
        sam=expand("mapped_reads/{sample}_{replicate}.sam",
            sample=SAMPLES, replicate=[1,2])
    output:
        "circexplorer2/back_spliced_junction.bed"
    log:
        "logs/circexplorer2/parse.log"
    params:
        aligner="BWA"
    # message:
    #     "CircExplorer2: extracting back-spliced exon-exon junction information"
    conda:
        "envs/circexplorer2.yaml"
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
    conda:
        "envs/circexplorer2.yaml"
    shell:
        "CIRCexplorer2 annotate -r {input.gene} -g {input.index} -b {input.bsj} "
        "-o {output} 2> {log}"
