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
# Last modification : 12-05-2021
################################################################################
import subprocess
import yaml

# VARIABLES
FAI_INDEX    = f'{PATH_genome}/{GENOME}.fna.fai'
HISAT2_INDEX = expand("{path}/hisat2/{genome}.{ext}.ht2", path = PATH_genome,
  genome = GENOME, ext = [1,2,3,4,5,6,7,8])

HISAT2_INDEX_QUANT = config["quantification"][["hisat2_index"]
BWA_INDEX_QUANT    = config["quantification"][["bwa_index"]
FAI_INDEX_QUANT    = config["quantification"]["fai_index"]
REFERENCE_QUANT    = config["quantification"]["reference"]
ANNOTATION_QUANT   = config["quantification"]["annotation"]



# TARGET RULE
rule quantification_results:
    input:
        circular = expand("{outdir}/ciriquant/{sample}.gtf", outdir = OUTDIR,
            sample = SAMPLES)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GENOME_INDEXS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule hisat2_index:
    output:
      HISAT2_INDEX
    params:
        prefix = GENOME,
        path   = PATH_genome
    priority: 10
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
    priority: 10
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
        config = f'{OUTDIR}/ciriquant/ciriquant_config.yaml'
    params:
        script = "src/utils/creating_yaml_file.py",
        genome = GENOME,
        hisat2 = True if HISAT2_INDEX_QUANT is not None else False,
        bwa    = True if BWA_INDEX_QUANT is not None else False,
        ref    = True if REFERENCE_QUANT is not None else False,
        gtf    = True if ANNOTATION_QUANT is not None else False
    conda: config["envs"]["ciriquant"]
    priority: 9
    shell:
        "python {params.script} {input.hisat2} {input.bwa} {input.ref} {input.gtf} {output}\
                {params.hisat2} {params.bwa} {params.fai} {params.ref} {params.gtf}"


rule ciriquant:
    input:
        read1  = lambda wildcards: f'{config["quantification"]["reads"]}/{wildcards.sample}{config["quantification"]["suffix"][1]}',
        read2  = lambda wildcards: f'{config["quantification"]["reads"]}/{wildcards.sample}{config["quantification"]["suffix"][2]}',
        ciri2  = lambda wildcards: f'{OUTDIR}/identification/overlap/{wildcards.sample}_common.txt', # potser aquest path també el podria facilitar l'usuari
        config = f'{OUTDIR}/ciriquant/ciriquant_config.yaml'
    output:
        # Linear transcripts alignment
        sorted_bam  = f'{OUTDIR}/ciriquant/align/{{sample}}.sorted.bam',
        sam_index   = f'{OUTDIR}/ciriquant/align/{{sample}}.sorted.bam.bai',
        # Linear gene abundance
        linear_transcripts_name = f'{OUTDIR}/ciriquant/gene/{{sample}}_cov.gtf', # StringTie outputs a file with the given name with all transcripts in the provided reference file that are fully covered by reads
        linear_out              = f'{OUTDIR}/ciriquant/gene/{{sample}}_out.gtf',
        gene_abundance          = f'{OUTDIR}/ciriquant/gene/{{sample}}_genes.list',
        # circRNAs quantification
        circular = f'{OUTDIR}/ciriquant/{{sample}}.gtf'
    threads: config["quantification"]["threads"]["ciriquant"]
    params:
        tool   = "CIRI2",
        outdir = f'{OUTDIR}/ciriquant',
    conda: config["envs"]["ciriquant"]
    priority: 8
    shell:
        "CIRIquant -t {threads} -1 {input.read1} -2 {input.read2} \
                    --config {input.config}\
                     -o {params.outdir}\
                     -p {wildcards.sample}\
                     --circ {input.ciri2}\
                     --tool {params.tool}"
