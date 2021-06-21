#!bin/bash/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "IN PROCESS"
################################################################################
# Snakefile to visualize improve the quantification of circRNA counts through
# CIRIquant tool.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              :
# Last modification : 01-06-2021
################################################################################
import subprocess
import yaml

# VARIABLES
## Generated
FAI_INDEX     = f'{PATH_genome}/{GENOME}.fna.fai'
HISAT2_INDEX  = expand("{path}/hisat2/{genome}.{ext}.ht2", path = PATH_genome,
  genome = GENOME, ext = [1,2,3,4,5,6,7,8])
PREFIX_HISAT2 = f'{PATH_genome}/hisat2/{GENOME}'
RES_QUANT     = f'{OUTDIR}/quantification/summary/summary.bed'
## Config file
CHECK_QUANT        = config["quantification"]["reads"]
CHECK_QUANT2       = config["quantification"]["circRNAs"]

R1_QUANT           = lambda wildcards: f'{config["quantification"]["reads"]}/{wildcards.sample}{config["quantification"]["read_suffix"][1]}'
R2_QUANT           = lambda wildcards: f'{config["quantification"]["reads"]}/{wildcards.sample}{config["quantification"]["read_suffix"][2]}'
CIRC_QUANT         = lambda wildcards: f'{config["quantification"]["circRNAs"]}/{wildcards.sample}{config["quantification"]["circ_suffix"]}'

HISAT2_INDEX_QUANT = config["quantification"]["hisat2_index"]
BWA_INDEX_QUANT    = config["quantification"]["bwa_index"]
FAI_INDEX_QUANT    = config["quantification"]["fai_index"]
REFERENCE_QUANT    = config["quantification"]["reference"]
ANNOTATION_QUANT   = config["quantification"]["annotation"]
CF_QUANT           = config["quantification"]["sample_threshold"]
CM_QUANT           = config["quantification"]["merged_threshold"]



# TARGET RULE
rule quantification_results:
    input:
        RES_QUANT

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GENOME_INDEXS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule hisat2_index:
    output:
      HISAT2_INDEX
    params:
        prefix = GENOME,
        path   = PATH_genome
    priority: 19
    shell:
      """
      genome={params.prefix}
      path={params.path}

      if [[ $genome == "GRCh38" ]];
        then
            echo "Downloading hisat2 index files from NCBI repository (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GO_TO_CURRENT_VERSION/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.hisat2_index.tar.gz)..."
            wget -c  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GO_TO_CURRENT_VERSION/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.hisat2_index.tar.gz \
                -O "$path/hisat2/GRCh38.hisat2_index.tar.gz"
            tar -zxvf "$path/hisat2/GRCh38.hisat2_index.tar.gz" -C "$path/hisat2/"

            pushd "$path/hisat2/"
            rm GRCh38.hisat2_index.tar.gz
            rename 's/GCA_000001405.15_GRCh38_no_alt_analysis_set/GRCh38/' *
            popd

      elif [[ $genome == "GRCh37" ]];
        then
            echo "Downloading hisat2 index files from UCSC repository (https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.hisat2_index.tar.gz)..."
            wget -c  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.hisat2_index.tar.gz \
                -O '$path/hisat2/GRCh37.hisat2_index.tar.gz'
            tar -zxvf "$path/hisat2/GRCh37.hisat2_index.tar.gz" -C "$path/hisat2/"

            pushd "$path/hisat2/"
            rm GRCh37.hisat2_index.tar.gz
            rename 's/hg19.p13.plusMT.no_alt_analysis_set/GRCh37/' *
            popd
        else
            echo "Please specify a genome in the config.yaml file. Genome: GRCh38(hg38) or GRCh37(hg37/hg19)."
        fi
      """

rule fai_index:
    # input:
    #     ref = f'{PATH_genome}/{GENOME}.fna'
    output:
        FAI_INDEX
    params:
        prefix = GENOME,
        path   = PATH_genome
    conda: config["envs"]["ciriquant"]
    priority: 19
    shell:
        """
        genome={params.prefix}
        path={params.path}

        if [[ $genome == "GRCh38" ]];
          then
            echo "Downloading a fai index from NCBI repository (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GO_TO_CURRENT_VERSION/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai)..."
            wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GO_TO_CURRENT_VERSION/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai\
                -O "$path/GRCh38.fna.fai"

        elif [[ $genome == "GRCh37" ]];
          then
            echo "Creating a fai index.."
            samtools faidx {input} -o {output}
        fi
        """

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CIRIQUANT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule ciriquant_config:
    input:
        hisat2 = HISAT2_INDEX_QUANT if HISAT2_INDEX_QUANT is not None else HISAT2_INDEX,
        bwa    = BWA_INDEX_QUANT if BWA_INDEX_QUANT is not None else BWA_INDEX,
        fai    = FAI_INDEX_QUANT if FAI_INDEX_QUANT is not None else FAI_INDEX,
        ref    = REFERENCE_QUANT if REFERENCE_QUANT is not None else REFERENCE,
        gtf    = ANNOTATION_QUANT if ANNOTATION_QUANT is not None else ANNOTATION
    output:
        config = f'{OUTDIR}/quantification/ciriquant_config.yaml'
    params:
        script        = "utils/creating_yaml_file.py",
        genome        = GENOME,
        prefix_bwa    = BWA_INDEX_QUANT if BWA_INDEX_QUANT is not None else PREFIX_BWA,
        prefix_hisat2 = HISAT2_INDEX_QUANT if HISAT2_INDEX_QUANT is not None else PREFIX_HISAT2
    conda: config["envs"]["ciriquant"]
    priority: 20
    shell:
        "python {params.script} {params.genome} {params.prefix_hisat2} {params.prefix_bwa} {input.ref} {input.gtf} {output}"


rule ciriquant:
    input:
        read1  = R1_QUANT if CHECK_QUANT is not None else TRI_R1 ,
        read2  = R2_QUANT if CHECK_QUANT is not None else TRI_R2,
        ciri2  = CIRC_QUANT if CHECK_QUANT2 is not None else RES_ID,
        config = f'{OUTDIR}/quantification/ciriquant_config.yaml'
    output:
        # Linear transcripts alignment
        sorted_bam  = f'{OUTDIR}/quantification/align/{{sample}}.sorted.bam',
        sam_index   = f'{OUTDIR}/quantification/align/{{sample}}.sorted.bam.bai',
        # Linear gene abundance
        linear_transcripts_name = f'{OUTDIR}/quantification/gene/{{sample}}_cov.gtf', # StringTie outputs a file with the given name with all transcripts in the provided reference file that are fully covered by reads
        linear_out              = f'{OUTDIR}/quantification/gene/{{sample}}_out.gtf',
        gene_abundance          = f'{OUTDIR}/quantification/gene/{{sample}}_genes.list',
        # circRNAs quantification
        circular = f'{OUTDIR}/quantification/{{sample}}.gtf'
    threads: config["quantification"]["threads"]["ciriquant"]
    params:
        tool    = "CIRI2",
        library = config["quantification"]["library_type"],
        outdir  = f'{OUTDIR}/quantification'
    conda: config["envs"]["ciriquant"]
    priority: 21
    shell:
        "CIRIquant -t {threads} -1 {input.read1} -2 {input.read2} \
                    --config {input.config}\
                     -o {params.outdir}\
                     -l {params.library}\
                     -p {wildcards.sample}\
                     --circ {input.ciri2}\
                     --tool {params.tool}"

rule ciriquant_results:
    input:
        expand("{outdir}/quantification/{sample}.gtf", outdir = OUTDIR, sample = SAMPLES)
    output:
        RES_QUANT
    params:
        tool   = "ciriquant",
        script = "utils/circM.py",
        sample_threshold = CF_QUANT,
        merged_threshold = CM_QUANT
    priority: 22
    shell:
        "python2 {params.script} -f {input} -a {params.tool}\
                -cf {params.sample_threshold} -cm {params.merged_threshold} > {output}"
