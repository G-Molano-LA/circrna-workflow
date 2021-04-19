#!/usr/bin/python3

from snakemake.utils import min_version
min_version("6.0")

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "ALMOST FINISHED" # requires execution to finish it

# 1. READING DATA###############################################################
data= [] # creating an empty list
with open("data/seqs.txt") as names:
    for line in names:
        data.append(line.rstrip('\n'))
SAMPLES=data
###############################################################################

rule results:
    input:
        results = expand("libs/identification/{sample}_coincident_circRNAs.txt",
            sample=SAMPLES),
        diagrams = expand("libs/identification/{sample}_venn_diagram.png",
            samples = SAMPLES)

"""
You need to have at least one rule (target rule) that does not produce any output,
but takes as input all your expected output. Snakemake will take that rule as
first rule, and then checks which rules produce the input of the rule. Then it
checks, which rules produce the inputs for those rules, etc.
"""
# 2. ALIGNMENT AND IDENFITICATION #############################################

rule dw_ref_genome:
    output:
        "data/raw_data/GRCh38.fa"
    message:
        " Downloanding reference genome (GRCh38)..."
    shell:
        "wget -c http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz "
        "-O data/raw_data/GRCh38.fa.gz && gunzip data/raw_data/GRCh38.fa.gz"

rule dw_ref_annotation:
    output:
        "data/raw_data/GRCh38_ann.gff"
    message:
        "Downloanding annotation reference ..."
    shell:
        "wget -c https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz "
        "-O data/raw_data/GRCh38_ann.gff.gz && gunzip data/raw_data/GRCh38_ann.gff.gz "

rule bwa_index:
    input:
        "data/raw_data/GRCh38.fa"
    output:
        expand("data/raw_data/bwa/{genome}.fna.{ext}", genome=["GRCh38"],
            ext=["amb", "ann", "bwt", "pac", "sa"])
    params:
        prefix="GRCh38",
        algorithm="bwtsw"
    message:
        "Creating a reference genome index from {params.prefix} file"
    log:
         "logs/bwa_index/GRCh38.log"
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa index -a {params.algorithm} -p {params.prefix} {input} 1> {output} 2> {log}"

rule mv_index:
    input:
        expand("{genome}.fa.{ext}", genome=["chr1"],
            ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        expand("data/raw_data/{genome}.fa.{ext}", genome=["chr1"],
            ext=["amb", "ann", "bwt", "pac", "sa"])
    shell:
        "mv {input} {output}"

rule bwa_mem:
    input:
        read1=lambda wildcards: "data/trimmed_data/"+wildcards.sample+"_1_val_1.fq.gz",
        read2=lambda wildcards: "data/trimmed_data/"+wildcards.sample+"_2_val_2.fq.gz"
    output:
        sam="data/mapped_data/{sample}.sam"
    params:
        score="19", # Do not output alignment with score lower than INT.
        prefix="data/raw_data/bwa/GRCh38"
    message:
        "Executing bwa_mem aligner with {threads} threads on the following files "
        "{input}. All alignment with score lower than {params.score} won't be output."
    threads: 10
    log:
        "logs/mapped_data/{sample}.log"
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa mem -T {params.score} -t {threads} {params.prefix} {input.read1} {input.read2} 1> {output.sam} 2> {log}"

rule dw_ciri2:
    output:
        "libs/identification/ciri2/CIRI2.pl"
    message:
        " Downloanding alignment tool CIRI2..."
    shell:
        ". src/utils/install_ciri2.sh"

rule ciri2_id:
    input:
        ciri="libs/identification/ciri2/CIRI2.pl",
        sam="data/mapped_data/{sample}.sam",
        ref="data/raw_data/GRCh38.fa" # no acepta archivo comprimido
    output:
        "libs/identification/ciri2/{sample}_results"
    message:
        "CIRI2:starting circRNA identification in {sample} file"
    log:
        "logs/ciri2/{sample}.log"
    conda:
        "envs/ciri2.yaml"
    shell:
        "perl {input.ciri} -I {input.sam} -O {output} -F {input.ref} 2> {log}"


rule circexplorer2_id:
    input:
        sam="data/mapped_data/{sample}.sam"
    output:
        "libs/identification/circexplorer2/{sample}_back_spliced_junction.bed"
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
        bsj="libs/identification/circexplorer2/{sample}_back_spliced_junction.bed",
        ref="data/raw_data/GRCh38.fa",
        gene="data/raw_data/GRCh38_ann.gff"
    output:
        "libs/identification/circexplorer2/{sample}_circularRNA_known.txt"
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

rule select_coincidences:
    input:
        ciri2="libs/identification/ciri2/{sample}_results",
        circexplorer2= "libs/identification/circexplorer2/{sample}_circularRNA_known.txt"
    output:
        "libs/identification/{sample}_coincident_circRNAs.txt",
        "libs/identification/{sample}_venn_diagram.png"
    shell:
        "Rscript src/utils/select_coincidents.R \
         --ciri2 {input.ciri2} \
         --circexplorer2 {input.circexplorer2}"

rule end:
    input:
        results = expand("libs/identification/{sample}_coincident_circRNAs.txt",
            sample=SAMPLES),
        diagrams = expand("libs/identification/{sample}_venn_diagram.png",
            samples = SAMPLES)
    output:
        "generated.txt"
    shell:
        "printf {input.results} > {output}"
