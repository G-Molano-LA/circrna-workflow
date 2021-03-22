#!bin/bash/python3

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
        "libs/ciriquant/ciri_quantification.gtf",
        "libs/ciriquant/circexp_quantification.gtf"

rule ciri2_quant:
    input:
        read1=expand("data/trimmed_data/{sample}_1_val_1.fq.gz",sample=SAMPLES),
        read2=expand("data/trimmed_data/{sample}_2_val_2.fq.gz",sample=SAMPLES)
    output:
        "libs/ciriquant/ciri_quantification.gtf"
    threads: 4
    params:
        tool="CIRI2",
        config="config.yaml", # falta configurar esto
        output_dir="libs/ciriquant/",
        pred_results="libs/ciri2/identification_results",
        prefix="ciri_quantification"
    conda:
        "envs/ciriquant.yaml"
    shell:
        "CIRIquant -t {threads} -1 {input.read1} -2 {input.read2} --config "
        " {params.config} -o {params.output_dir} -p {params.prefix} "
        " --circ {params.pred_results} --tool {params.tool}"

rule circexp_quant:
    input:
        read1=expand("data/trimmed_data/{sample}_1_val_1.fq.gz",sample=SAMPLES),
        read2=expand("data/trimmed_data/{sample}_2_val_2.fq.gz",sample=SAMPLES)
    output:
        "libs/ciriquant/circexp_quantification.gtf"
    threads: 4
    params:
        tool="CIRCexplorer2",
        config="config.yaml", # falta configurar esto
        output_dir="libs/ciriquant/",
        pred_results="libs/circexplorer2/circularRNA_known.txt",
        prefix="circexp_quantification"
    conda:
        "envs/ciriquant.yaml"
    shell:
        "CIRIquant -t {threads} -1 {input.read1} -2 {input.read2} --config "
        " {params.config} -o {params.output_dir} -p {input.prefix} "
        " --circ {params.pred_results} --tool {params.tool}"
