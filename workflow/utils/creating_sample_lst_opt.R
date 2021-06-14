#!/bin/R

################################################################################
# Rscript to obtain configuration files (containing sample name, sample path and
# condition samples such as group and sex) that are required for CIRIquant_prep.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 04-05-2021
# Last modification : 10-05-2021
################################################################################

suppressPackageStartupMessages(library("argparse"))

source("utils/utils.R")

# SCRIPT ARGUMENTS
parser <- ArgumentParser(description = 'Create some list')
parser$add_argument("--metadata", required = TRUE,  default = NULL, action = 'store',
  help = "File containing sample names and metadata variables such as group and sex.")
parser$add_argument("--sep", default = ",", action = 'store', help = "Separator")
parser$add_argument("--dir",  required = TRUE, action = 'store', default = NULL,
  help = "directory")
parser$add_argument("--outdir", action = 'store', default = NULL, help = "output directory")

opt <- parser$parse_args()

if(is.null(opt$metadata) || is.null(opt$sep) || is.null(opt$dir)){
  parser$print_help()
  stop("Options --metadata/--sep/--dir must be supplied\n", call. = FALSE)
}


# Passing args
sep      <- check_sep(opt$sep)
metadata <- check_metadata(opt$metadata, sep)
dir      <- opt$dir
outdir   <- opt$outdir


# Creating content
# 1. sample.lst
circular_path <- vector()
for (i in 1:length(metadata$sample)){
  circular_path[i] <- paste0(dir,"/",metadata$sample[i],".gtf")}

if('sex' %in% colnames(metadata)){
  sample_lst <- cbind(metadata$sample, circular_path, as.vector(metadata$group), as.vector(metadata$sex))
}else{
  sample_lst <- cbind(metadata$sample, circular_path, as.vector(metadata$group))
}


# 2. sample_gen.lst
linear_dir <- paste0(dir,"/gene/")

linear_path <- vector()
for (i in 1:length(metadata$sample)){
  linear_path[i] <- paste0(linear_dir,metadata$sample[i],"_out.gtf")}

sample_gene_lst <- cbind(metadata$sample, linear_path)


# Write tsv file
write.table(sample_lst, file = paste0(outdir,"/sample_circ.lst"), row.names = FALSE,
  col.names = FALSE, quote = FALSE)
write.table(sample_gene_lst, file = paste0(outdir,"/sample_gene.lst"), row.names = FALSE,
  col.names = FALSE, quote = FALSE)
