author = "G. Molano, LA (gonmola@hotmail.es)"

# Passing arguments
args = commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
}else{
  outdir = args[1]
}


# Open metadata file
metadata <- read.csv("libs/DE_analysis/metadata_samples.csv")

# Selecting desired information
sample_names <- metadata$Run

# Creating content
# 1. sample.lst
circular_path <- vector()
for (i in 1:length(sample_names)){
  circular_path[i] <- paste0(outdir,"/",sample_names[i],".gtf")}

sample_lst <- cbind(sample_names, circular_path, metadata["group"], metadata["sex"])

# 2. sample_gen.lst
linear_dir <- paste0(outdir,"/gene/")

linear_path <- vector()
for (i in 1:length(sample_names)){
  linear_path[i] <- paste0(linear_dir,sample_names[i],"_out.gtf")}

sample_gene_lst <- cbind(sample_names, linear_path)


# Write tsv file
write.table(sample_lst, file = "libs/ciriquant/sample.lst", row.names = FALSE,
  col.names = FALSE, quote = FALSE)
write.table(sample_gene_lst, file = "libs/ciriquant/sample_gene.lst", row.names = FALSE,
  col.names = FALSE, quote = FALSE)
