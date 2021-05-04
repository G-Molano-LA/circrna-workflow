#!/bin/R

suppressPackageStartupMessages(library("argparse"))
author = "G. Molano, LA (gonmola@hotmail.es)"


parser <- ArgumentParser(description = 'Create some list')
parser$add_argument("--samples", action = 'store', nargs = '+', default = NULL,
  help = "list of sample names")
parser$add_argument("--group", action = 'store', nargs = '+', default = NULL,
  help = "list of group values")
parser$add_argument("--sex", action = 'store', nargs = '+', default = NULL,
  help = "list of sex values")
parser$add_argument("--dir", action = 'store', default = NULL, help = "directory")
parser$add_argument("--outdir", action = 'store', default = NULL, help = "output directory")

opt <- parser$parse_args()

if(is.null(opt$samples) || is.null(opt$dir) || is.null(opt$group)){
  parser$print_help()
  stop("Options --samples/--dir/--group must be supplied\n", call. = FALSE)
}


# Passing args
sample_names <- opt$samples
sample_group <- opt$group
sample_sex   <- opt$sex
dir          <- opt$dir
outdir       <- opt$outdir


# Creating content
# 1. sample.lst
circular_path <- vector()
for (i in 1:length(sample_names)){
  circular_path[i] <- paste0(dir,"/",sample_names[i],".gtf")}

sample_lst <- cbind(sample_names, circular_path, sample_group, sample_sex)

# 2. sample_gen.lst
linear_dir <- paste0(dir,"/gene/")

linear_path <- vector()
for (i in 1:length(sample_names)){
  linear_path[i] <- paste0(linear_dir,sample_names[i],"_out.gtf")}

sample_gene_lst <- cbind(sample_names, linear_path)


# Write tsv file
write.table(sample_lst, file = paste0(outdir,"/DE_analysis/prep_DE/sample_circ.lst"), row.names = FALSE,
  col.names = FALSE, quote = FALSE)
write.table(sample_gene_lst, file = paste0(outdir,"/DE_analysis/prep_DE/sample_gene.lst"), row.names = FALSE,
  col.names = FALSE, quote = FALSE)
