#!/usr/bin/python3

__author__ = "G. Molano, LA (gonmola@hotmail.es)"
__state__ = "ALMOST FINISHED" # Falta comprovar funcionamiento con GRCh37

################################################################################
# Snakefile to align and identification circRNAs from RNA-seq data.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              :
# Last modification : 01-06-2021
################################################################################

# VARIABLES
## Generated
REFERENCE    = f'{PATH_genome}/{GENOME}.fna'
ANNOTATION   = f'{PATH_genome}/{GENOME}_ann.gtf'
REFFLAT_ANN  = f'{PATH_genome}/{GENOME}.refFlat.txt'
BWA_INDEX    = expand("{path}/bwa/{genome}.fna.{ext}", path = PATH_genome,
    genome = GENOME, ext=["amb", "ann", "bwt", "pac", "sa"])
PREFIX_BWA   = f'{PATH_genome}/bwa/{GENOME}.fna'
RES_ID       = lambda wildcards: f'{OUTDIR}/identification/overlap/{wildcards.sample}_common.txt'
## Config file
CHECK_ALN   = config["aln_and_id"]["reads"]

BWA_ALN     = config["aln_and_id"]["bwa_index"]
REF_ALN     = config["aln_and_id"]["reference"]
ANN_ALN     = config["aln_and_id"]["annotation"]
REFFLAT_ALN = config["aln_and_id"]["refFlat"]
CF          = config["aln_and_id"]["sample_threshold"]
CM          = config["aln_and_id"]["merged_threshold"]
R1_ALN      = lambda wildcards: f'{config["aln_and_id"]["reads"]}/{wildcards.sample}{config["aln_and_id"]["suffix"][1]}'
R2_ALN      = lambda wildcards: f'{config["aln_and_id"]["reads"]}/{wildcards.sample}{config["aln_and_id"]["suffix"][2]}'


# TARGET RULE
rule alignment_and_identification_results:
    input:
        bed  = f'{OUTDIR}/identification/overlap/summary_overlap.bed',
        venn = f'{OUTDIR}/identification/overlap/summary_overlap.png'

#~~~~~~~~~~~~~~~~~~~~~~~~~ANNOTATION&REFERENCE-FILES~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule get_ref_genome:
    output:
        ref = REFERENCE
    params:
        path = PATH_genome
    priority: 95
    # message:
    #     "Downloanding"
    run:
        if GENOME == 'GRCh38':
            shell('wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GO_TO_CURRENT_VERSION/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz \
                            -O {params.path}/GRCh38.fna.gz \
                            && gunzip {params.path}/GRCh38.fna.gz',
                            shell = True, text = True) # UCSC format
        elif GENOME == 'GRCh37':
            shell('wget -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz \
                            -O {params.path}/GRCh37.fna.gz \
                            && gunzip {params.path}/GRCh37.fna.gz',
                            shell = True, text = True) # ensEMBL format
        else:
            print("Please specify a genome in the config.yaml file. Genome: hg38/GRCh38 or hg37/GRCh37/hg19")


rule get_ref_annotation:
    output:
        ANNOTATION
    params:
        path = PATH_genome
    priority: 94
    run:
        if GENOME == 'GRCh38':
            shell("wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz \
                            -O {params.path}/GRCh38_ann.gtf.gz \
                            && gunzip {params.path}/GRCh38_ann.gtf.gz",
                            shell = True, text = True) # UCSC format
        elif GENOME== 'GRCh37':
            shell("wget -c http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz \
                            -O {params.path}/GRCh37_ann.gtf.gz \
                            && gunzip {params.path}/GRCh37_ann.gtf.gz",
                            shell = True, text = True) # ensEMBL format


rule refFlat:
    output:
        REFFLAT_ANN
    params:
        path = PATH_genome
    priority: 93
    run:
        if GENOME == 'GRCh38':
            shell("wget -c https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refFlat.txt.gz \
                            -O {params.path}/GRCh38_refFlat.txt.gz \
                            && gunzip {params.path}/GRCh38_refFlat.txt.gz",
                            shell = True, text = True)
        elif GENOME == 'GRCh37':
            shell("bash src/utils/gtf_to_refFlat.sh {params.path}",
                            shell = True, text = True)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~BWA-MEM ALIGNMENT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule bwa_index:
    input:
        REFERENCE
    output:
        BWA_INDEX
    params:
        genome = GENOME,
        script = "src/utils/bwa_index.sh",
        path = PATH_genome
    conda: config["envs"]["bwa"]
    priority: 93
    shell:
        """
        genome={params.genome}
        path={params.path}

        if [[ $genome == "GRCh38" ]];
        then
            echo "Downloading genome index files from NCBI repository (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GO_TO_CURRENT_VERSION/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz)..."
            wget -c  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GO_TO_CURRENT_VERSION/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz \
                -O "$path/bwa/GRCh38.bwa_index.tar.gz"
            tar -zxvf "$path/bwa/GRCh38.bwa_index.tar.gz" -C "$path/bwa"

            pushd "$path/bwa/"
            rm GRCh38.bwa_index.tar.gz
            rename 's/GCA_000001405.15_GRCh38_no_alt_analysis_set/GRCh38/' *
            popd
        elif [[ $genome == "GRCh37" ]];
        then
            echo "Creating a genome index file from GRCh37.fna genome"
            bwa index -a bwtsw -p "$genome" "$path/$genome.fna"
        else
            echo "Please specify a genome in the config.yaml file. Genome: hg38/GRCh38 or hg37/GRCh37/hg19."
        fi
        """

rule bwa_mem:
    input:
        read1 = R1_ALN if CHECK_ALN is not None else TRI_R1,
        read2 = R2_ALN if CHECK_ALN is not None else TRI_R2,
        index = BWA_ALN if BWA_ALN is not None else BWA_INDEX
    output:
        sam    = temp(f'{OUTDIR}/data/mapped_data/{{sample}}.sam')
    params:
        score  = "19", # Do not output alignment with score lower than INT.
        prefix = BWA_ALN if BWA_INDEX is not None else PREFIX_BWA
    message:
        "BWA-MEM aligner: Starting the alignment of the reads from {input.read1} & {input.read2}\
         THREADS = {threads}\
         OUTPUT  = {output}\
         GENOME  = {params.prefix}\
         SCORE   = All alignment with score lower than {params.score} won't be output\
         LOG     = {log}"
    threads: config["aln_and_id"]["threads"]["bwa"]
    log:
        f'{OUTDIR}logs/mapped_data/{{sample}}.log'
    conda: config["envs"]["bwa"]
    priority: 92
    shell:
        "bwa mem -T {params.score} -t {threads} {params.prefix} {input.read1} {input.read2} 1> {output.sam} 2> {log}"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CIRI2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule dw_ciri2:
    output:
        "src/tools/CIRI2.pl"
    message:
        "Downloanding alignment tool CIRI2 from official website...\
        SCRIPT = {params.script}"
    params:
        script = "src/utils/install_ciri2.sh"
    priority: 91
    shell:
        """
        wget https://sourceforge.net/projects/ciri/files/CIRI2/CIRI_v2.0.6.zip
        unzip CIRI_v2.0.6.zip
        rm -r CIRI_v2.0.6.zip
        mv CIRI_v2.0.6/CIRI2.pl src/tools/
        rm -r CIRI_v2.0.6/
        rm -r __MACOSX/
        """


rule ciri2:
    input:
        ciri = "src/tools/CIRI2.pl",
        sam  = f'{OUTDIR}/data/mapped_data/{{sample}}.sam',
        ref  =  REF_ALN if REF_ALN is not None else REFERENCE # no acepta archivo comprimido
    output:
        f'{OUTDIR}/identification/ciri2/{{sample}}_results.txt'
    threads: config["aln_and_id"]["threads"]["ciri2"]
    message:
        "CIRI2: Starting circRNA identification in {wildcards.sample}.sam file...\
        REFERENCE FILE = {input.ref}\
        THREADS        = {threads} \
        OUTPUT         = {output}"
    conda: config["envs"]["ciri2"]
    priority: 90
    shell:
        "perl {input.ciri} -I {input.sam} -O {output} -F {input.ref} -T {threads}"

rule ciri2_results:
    input:
        expand("{outdir}/identification/ciri2/{sample}_results.txt",
            outdir = OUTDIR, sample = SAMPLES)
    output:
        f'{OUTDIR}/identification/ciri2/ciri2_results.bed'
    params:
        tool             = "ciri2",
        script           = "src/utils/circM.py",
        sample_threshold = CF,
        merged_threshold = CM
    priority: 89
    conda: config["envs"]["R"]
    shell:
        "python2 {params.script} -f {input} -a {params.tool}\
                -cf {params.sample_threshold} -cm {params.merged_threshold} > {output}"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CIRCEXPLORER2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule circexplorer2_id:
    input:
        sam = f'{OUTDIR}/data/mapped_data/{{sample}}.sam',
    output:
        temp(f'{OUTDIR}/identification/circexplorer2/{{sample}}_back_spliced_junction.bed')
    params:
        aligner="BWA"
    message:
        "CIRCEXPLORER2: Extracting back-spliced exon-exon junction information from INPUT = {input} ...\
        ALIGNER = {params.aligner}\
        OUTPUT  = {output}"
    conda: config["envs"]["circexplorer2"]
    priority: 90
    shell:
        "CIRCexplorer2 parse -t {params.aligner} --bed={output} {input}"


rule circexplorer2_annotation:
    input:
        bsj     = f'{OUTDIR}/identification/circexplorer2/{{sample}}_back_spliced_junction.bed',
        ref     = REF_ALN if REF_ALN is not None else REFERENCE,
        refFlat = REFFLAT_ALN if REFFLAT_ALN is not None else REFFLAT_ANN
    output:
        f'{OUTDIR}/identification/circexplorer2/{{sample}}_circularRNA_known.txt'
    message:
        "CIRCEXPLORER2: Annotating circRNAs with known RefSeq genes in {input.bsj}. \
        REFERENCE FILE  = {input.ref}\
        ANNOTATING FILE = {input.refFlat} (refFlat format) \
        OUTPUT          = {output}"
    conda: config["envs"]["circexplorer2"]
    priority: 89
    shell:
        "CIRCexplorer2 annotate -r {input.refFlat} -g {input.ref} -b {input.bsj} "
        "-o {output}"

rule circexplorer2_results:
    input:
        expand("{outdir}/identification/circexplorer2/{sample}_circularRNA_known.txt",
            outdir = OUTDIR, sample = SAMPLES)
    output:
        f'{OUTDIR}/identification/circexplorer2/circexplorer2_results.bed'
    params:
        tool   = "circexplorer2",
        script = "src/utils/circM.py",
        sample_threshold = CF,
        merged_threshold = CM
    priority: 88
    shell:
        "python2 {params.script} -f {input} -a {params.tool}\
                -cf {params.sample_threshold} -cm {params.merged_threshold} > {output}"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OVERLAP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rule select_coincidences:
    input:
        ciri2         = f'{OUTDIR}/identification/ciri2/{{sample}}_results.txt',
        circexplorer2 = f'{OUTDIR}/identification/circexplorer2/{{sample}}_circularRNA_known.txt'
    output:
        txt  = f'{OUTDIR}/identification/overlap/{{sample}}_common.txt',
        venn = f'{OUTDIR}/identification/overlap/{{sample}}_common.png'
    params:
        script = "src/utils/select_coincidents.R",
        outdir = f'{OUTDIR}/identification/overlap/',
        sample = f'{{sample}}'
    message:
        "OVERLAP: Selecting circRNAs commonly identified by CIRI2 and CircExplorer2 tools... \
        INPUTS  = 1) {input.ciri2} ;  2){input.circexplorer2} \
        OUTPUTS = 1) Venn Diagram: {output.venn} ; 2) Circular matrix: {output.txt} \
        SCRIPT  = {params.script}"
    priority: 87
    conda: config["envs"]["R"]
    shell:
        "Rscript {params.script} \
         --ciri2 {input.ciri2} \
         --circexplorer2 {input.circexplorer2}\
         --sample {params.sample}\
         --outdir {params.outdir}"

rule overlap_results:
    input:
        expand("{outdir}/identification/overlap/{sample}_common.txt",
            outdir = OUTDIR, sample = SAMPLES)
    output:
        f'{OUTDIR}/identification/overlap/summary_overlap.bed'
    params:
        tool   = "ciri2",
        script = "src/utils/circM.py",
        sample_threshold = CF,
        merged_threshold = CM
    priority: 86
    shell:
        "python2 {params.script} -f {input} -a {params.tool}\
                -cf {params.sample_threshold} -cm {params.merged_threshold} > {output}"

rule overlap_visualization:
    input:
        ciri2         = f'{OUTDIR}/identification/ciri2/ciri2_results.bed',
        circexplorer2 = f'{OUTDIR}/identification/circexplorer2/circexplorer2_results.bed',
        overlap       = f'{OUTDIR}/identification/overlap/summary_overlap.bed'
    output:
        venn = f'{OUTDIR}/identification/overlap/summary_overlap.png'
    params:
        script = "src/utils/overlap_visualization.R",
        outdir = f'{OUTDIR}/identification/overlap/'
    message:
        "OVERLAP: Selecting circRNAs commonly identified by CIRI2 and CircExplorer2 tools... \
        INPUTS  = 1) {input.ciri2} ;  2){input.circexplorer2}; 3) {input.overlap} \
        OUTPUTS = 1) Venn Diagram: {output.venn} \
        SCRIPT  = {params.script}"
    priority: 85
    conda: config["envs"]["R"]
    shell:
        "Rscript {params.script} \
         --ciri2 {input.ciri2} \
         --circexplorer2 {input.circexplorer2}\
         --overlap {input.overlap}\
         --outdir {params.outdir}"
