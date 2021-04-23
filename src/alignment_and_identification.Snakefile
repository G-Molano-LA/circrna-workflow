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
        diagrams = expand("libs/plots/venn_diagrams/{sample}.png", sample = SAMPLES)

"""
You need to have at least one rule (target rule) that does not produce any output,
but takes as input all your expected output. Snakemake will take that rule as
first rule, and then checks which rules produce the input of the rule. Then it
checks, which rules produce the inputs for those rules, etc.
"""
# 2. ALIGNMENT AND IDENFITICATION #############################################
#
# rule bwa_ref_genome:
#     output:
#         "data/raw_data/GRCh38.fna"
#     priority: 10
#     shell:
#         "wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GO_TO_CURRENT_VERSION/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz \
#           -O data/raw_data/GRCh38.fna.gz \
#           && gunzip data/raw_data/GRCh38.fna.gz \
#           && rm data/raw_data/GRCh38.fna.gz"
#
# rule bwa_dw_index:
#     output:
#         expand("{genome}.fna.{ext}", genome=["GRCh38"], ext=["amb", "ann", "bwt", "pac", "sa"])
#     priority: 10
#     shell:
#         "wget -c  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.bwa_index.tar.gz \
#             -O data/raw_data/bwa/GCA_000001405.15_GRCh38_full_analysis_set.fna.bwa_index.tar.gz \
#           && tar -zxvf data/raw_data/bwa/GCA_000001405.15_GRCh38_full_analysis_set.fna.bwa_index.tar.gz -C data/raw_data/bwa \
#           && pushd data/raw_data/bwa/ \
#           && rm GCA_000001405.15_GRCh38_full_analysis_set.fna.bwa_index.tar.gz \
#           && rename 's/GCA_000001405.15_GRCh38_full_analysis_set/GRCh38/' * \
#           && popd"
#
# rule dw_ref_annotation:
#     output:
#         "data/raw_data/GRCh38_ann.gtf"
#     priority: 10
#     shell:
#         "wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz \
#         -O data/raw_data/GRCh38_ann.gtf.gz && gunzip data/raw_data/GRCh38_ann.gtf.gz"
#
# rule gtf_to_refFlat:
#     input:
#         "data/raw_data/GRCh38_ann.gtf"
#     output:
#         "data/raw_data/GRCh38_refFlat.txt"
#     params:
#         script = "src/utils/gtf_to_refFlat.sh"
#     conda: "envs/circexplorer2.yaml"
#     shell:
#         "bash {params.script}"

rule bwa_mem:
    input:
        read1=lambda wildcards: "data/trimmed_data/"+wildcards.sample+"_1_val_1.fq.gz",
        read2=lambda wildcards: "data/trimmed_data/"+wildcards.sample+"_2_val_2.fq.gz",
        index=expand("data/raw_data/bwa/{genome}.fna.{ext}", genome=["GRCh38"],
            ext=["amb", "ann", "bwt", "pac", "sa"])
    output:
        sam=temp("data/mapped_data/{sample}.sam")
    params:
        score="19", # Do not output alignment with score lower than INT.
        prefix="data/raw_data/bwa/GRCh38.fna"
    # message:
    #     "Executing bwa_mem aligner with {threads} threads on the following files "
    #     "{input}. All alignment with score lower than {params.score} won't be output."
    threads: 5
    log:
        "logs/mapped_data/{sample}.log"
    conda:
        "envs/bwa.yaml"
    priority: 7
    shell:
        "bwa mem -T {params.score} -t {threads} {params.prefix} {input.read1} {input.read2} 1> {output.sam} 2> {log}"

rule dw_ciri2:
    output:
        "libs/identification/ciri2/CIRI2.pl"
    # message:
    #     " Downloanding alignment tool CIRI2..."
    priority: 6
    shell:
        ". src/utils/install_ciri2.sh"

rule ciri2_id:
    input:
        ciri="libs/identification/ciri2/CIRI2.pl",
        sam="data/mapped_data/{sample}.sam",
        gene="data/raw_data/GRCh38_ann.gtf",
        ref="data/raw_data/GRCh38.fna" # no acepta archivo comprimido
                                        # Mirar qué tipo utiliza
    output:
        "libs/identification/ciri2/{sample}_results"
    # message:
    #     "CIRI2:starting circRNA identification in {sample} file"
    log:
        "logs/ciri2/{sample}.log"
    conda:
        "envs/ciri2.yaml"
    priority: 5
    shell:
        "perl {input.ciri} -I {input.sam} -O {output} -F {input.ref} - A {input.gene} 2> {log}"


rule circexplorer2_id:
    input:
        sam="data/mapped_data/{sample}.sam"
    output:
        "libs/identification/circexplorer2/{sample}_back_spliced_junction.bed"
    log:
        "logs/circexplorer2/{sample}_parse.log"
    params:
        aligner="BWA"
    # message:
    #     "CircExplorer2: extracting back-spliced exon-exon junction information from {input.sam}"
    conda:
        "envs/circexplorer2.yaml"
    priority: 6
    shell:
        "CIRCexplorer2 parse -t {params.aligner} --bed={output} {input} 2> {log}"

rule circexplorer2_annotation:
    input:
        bsj="libs/identification/circexplorer2/{sample}_back_spliced_junction.bed",
        ref="data/raw_data/GRCh38.fna",
        gene="data/raw_data/hg38_refFlat.txt"
    output:
        "libs/identification/circexplorer2/{sample}_circularRNA_known.txt"
    log:
        "logs/circexplorer2/{sample}_annotate.log"
    # message:
    #     "CircExplorer2: annotating circRNAs with known RefSeq genes in {input.bsj}. \
    #     The output file {output} will be generated."
    conda:
        "envs/circexplorer2.yaml"
    priority: 5
    shell:
        "CIRCexplorer2 annotate -r {input.gene} -g {input.ref} -b {input.bsj} "
        "-o {output} 2> {log}"

rule select_coincidences:
    input:
        ciri2="libs/identification/ciri2/{sample}_results",
        circexplorer2= "libs/identification/circexplorer2/{sample}_circularRNA_known.txt"
    output:
        "libs/identification/{sample}_coincident_circRNAs.txt",
        "libs/plots/venn_diagrams/{sample}.png"
    priority: 4
    conda :
        "envs/R.yaml"
    shell:
        "Rscript src/utils/select_coincidents.R \
         --ciri2 {input.ciri2} \
         --circexplorer2 {input.circexplorer2}"

rule end:
    input:
        results = expand("libs/identification/{sample}_coincident_circRNAs.txt",
            sample=SAMPLES),
        diagrams = expand("libs/plots/venn_diagrams/{sample}.png",
            sample = SAMPLES)
    output:
        "generated.txt"
    shell:
        "printf {input.results} > {output}"
