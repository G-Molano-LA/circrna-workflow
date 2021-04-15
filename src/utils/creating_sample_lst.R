author = "G. Molano, LA (gonmola@hotmail.es)"

# Open metadata file
metadata <- read.csv("libs/DE_analysis/metadata_samples.csv")

# Selecting desired information
sample_names <- metadata$Run

# Creating content
# 1. sample.lst
circ_dir <- "/libs/ciriquant/" # Revisar este directorio cuando se hayan seleccionado los coincidentes

circular_path <- vector()

for (i in 1:length(sample_names)){
  circular_path[i] <- paste0(circ_dir,sample_names[i],".gtf")}

sample_lst <- cbind(sample_names, circular_path, metadata["group"], metadata["sex"])

# 2. sample_gen.lst
linear_dir <- "/libs/ciriquant/gene/"

linear_path <- vector()

for (i in 1:length(sample_names)){
  linear_path[i] <- paste0(linear_dir,sample_names[i],"_out.gtf")}

sample_gene_lst <- cbind(sample_names, linear_path)


# Write tsv file
write.table(sample_lst, file = "libs/ciriquant/sample.lst", row.names = FALSE,
  col.names = FALSE)
write.table(sample_gene_lst, file = "libs/ciriquant/sample_gene.lst", row.names = FALSE,
  col.names = FALSE)
