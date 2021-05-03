

suppressPackageStartupMessages(library("optparse"))
author = "G. Molano, LA (gonmola@hotmail.es)"


option_list <- list(
  make_option(c("--samples"), default = NULL, type = "character", help = "sample names"),
  make_option(c("--dir"), default = NULL, type = "character", help = "sample directory"),
  make_option(c("--group"), default = NULL, type = "character", help "group variable. e.g control or treatment"),
  make_option(c("--sex"), default = NULL,type = "character", help = "sex variable"),
  make_option(c("--outdir"), default = NULL)
)

parser <- OptionParser(option_list = option_list)
opt    <- parser_args(parser)

if(is.null(opt$samples) || is.null(opt$dir) || is.null(opt$group)){
  print_help(parser)
  stop("Options --samples/--dir/--group must be supplied\n", call. = FALSE)
}

# Passing args

sample_names <- as.vector(opt$samples)
sample_group <- as.vector(opt$group)
sample_sex   <- as.vector(opt$sex)
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
