#!bin/bash/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "IN PROCESS"

# Depend on a minimum snakemake version
from snakemake.utils import min_version
min_version("6.0")
###############################################################################
data= [] # creating an empty list
with open("data/seqs.txt") as fp:
    for line in fp:
        data.append(line.rstrip('_[1,2].fastq.gz\n'))
data=list(dict.fromkeys(data)) # removing duplicates
SAMPLES=data
###############################################################################

rule results:
    input:
        "generated.txt"

rule hisat2_index:
    input:
       "data/raw_data/GRCh38.fa"
    output:
      expand("data/raw_data/hitsa2/{genome}.{ext}.ht2", genome=["GRCh38"],
            ext=[1,2,3,4,5,6,7,8])
    params:
        prefix = "GRCh38"
    threads: 16
    conda:
       "envs/ciriquant.yaml"
    shell:
      "hisat2-build -p {threads} {input} {params.prefix}"

rule reference_fai_index:
    input:
        "data/raw_data/GRCh38.fa"
    output:
        "data/raw_data/hisat2/GRCh38.fa.fai"
    conda:
        "envs/ciriquant.yaml"
    shell:
        "samtools faidx {input} -o {output}"

rule write_config_yaml:
    input:
        expand("data/raw_data/hitsa2/{genome}.{ext}.ht2", genome=["GRCh38"],
              ext=[1,2,3,4,5,6,7,8]),
        "data/raw_data/hisat2/GRCh38.fa.fai"
    output:
        config="libs/ciriquant/config_file.yaml"
    script:
        "src/utils/creating_yaml_file.py"

rule ciriquant:
    input:
        read1=lambda wildcards: "data/trimmed_data/"+wildcards.sample+"_1_val_1.fq.gz",
        read2=lambda wildcards: "data/trimmed_data/"+wildcards.sample+"_2_val_2.fq.gz",
        config="libs/ciriquant/config_file.yaml"
    output:
        circular="libs/ciriquant/ciri2/{sample}.gtf", # REVISAR si este output se genera en este directorio
        bam=temp("libs/ciriquant/align/{sample}.bam"),
        sorted_bam=temp("libs/ciriquant/align/{sample}.sorted.bam"),
        sam_index=temp("libs/ciriquant/align/{sample}.sorted.bam.bai"),
        linear_transcripts_name="libs/ciriquant/gene/{sample}_cov.gtf",
        linear="libs/ciriquant/gene/{sample}_out.gtf",
        gene_abundance="libs/ciriquant/gene/{sample}_genes.list"
    threads: 4
    params:
        tool="CIRI2",
        output_dir="libs/ciriquant/",
        pred_results = lambda wildcards : "libs/ciri2/"+wildcards.sample+"_results"
    conda:
        "envs/ciriquant.yaml"
    shell:
        "CIRIquant -t {threads} -1 {input.read1} -2 {input.read2} --config "
        " {input.config} -o {params.output_dir}"
        " --circ {params.pred_results} --tool {params.tool}"

rule end:
    input:
        ciri=expand("libs/ciriquant/ciri2/{sample}.gtf", sample=SAMPLES)
    output:
        "generated.txt"
    shell:
        "printf {input.ciri} > {output}"
