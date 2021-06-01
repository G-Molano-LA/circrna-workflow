from snakemake.utils import min_version
import subprocess

min_version("6.0")


# CONFIG FILE
configfile: "src/utils/config.yaml"

# VARIABLES
SAMPLES      = config["samples"]
PATH_genome  = config["genome_files"]
GENOME       = config["genome"]
OUTDIR       = config["outdir"]

# MODULES
QUALITY         = config['modules']['quality_control']
TRIMMING        = config['modules']['trimming']
ALN_and_ID      = config['modules']['alignment_and_identification']
ANNOTATION      = config["modules"]['annotation']
QUANTIFICATION  = config['modules']['quantification']
DE_ANALYSIS     = config['modules']['DE_analysis']
VISUALIZATION   = config['modules']['visualization']

# RULES
include: "src/rules/quality_control.smk"
include: "src/rules/trimming.smk"
include: "src/rules/alignment_and_identification.smk"
include: "src/rules/quantification.smk"
include: "src/rules/annotation.smk"
include: "src/rules/DE_analysis.smk"
include: "src/rules/visualization.smk"

rule all:
    input:
        # Initial quality control
        *(rules.quality_control_results.input if QUALITY else []),

        # trimming
        *(rules.trimming_results.input if TRIMMING  else []),

        # alignment_and_identification
        *(rules.alignment_and_identification_results.input if ALN_and_ID  else []),

        # quantification
        *(rules.quantification_results.input if QUANTIFICATION else []),

        # Annotation
        *(rules.annotation_results.input if ANNOTATION else []),

        # Differential expression analysis,
        *(rules.DE_results.input if QUANTIFICATION and DE_ANALYSIS else []),

        # Visualization
        *(rules.visualization_results.input if VISUALIZATION else [])
