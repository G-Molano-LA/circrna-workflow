from snakemake.utils import min_version
import subprocess

min_version("6.0")


# CONFIG FILE
configfile: "src/utils/config.yaml"

# VARIABLES
SAMPLES      = config["samples"]
PATH_genome  = config["path"]["genome_files"]
GENOME       = config["genome"]
OUTDIR       = config["path"]["outdir"]

# RULES
include: "src/rules/quality_control.smk"
include: "src/rules/trimming.smk"
include: "src/rules/alignment_and_identification.smk"
include: "src/rules/quantification.smk"
include: "src/rules/DE_analysis.smk"
include: "src/rules/visualization.smk"

rule all:
    input:
        # Initial quality control
        *(rules.quality_control_results.input if config['modules']['quality_control'] else []),

        # trimming
        *(rules.trimming_results.input if config['modules']['trimming'] else []),

        # alignment_and_identification
        *(rules.alignment_and_identification_results.input if config['modules']['alignment_&_identification'] else []),

        # quantification
        *(rules.quantification_results.input if config['modules']['quantification'] else []),

        # Differential expression analysis
        *(rules.quantification_results.input if config['modules']['quantification'] else []) # Doble condition (also if 'DE' yes)
