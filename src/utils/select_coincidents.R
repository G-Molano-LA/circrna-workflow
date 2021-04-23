#!/bin/R

###############################################################################
## R script for selecting common identified circRNAS by CIRI2 & CircExplorer2
## Author: G. Molano, LA (gonmola@hotmail.es)
###############################################################################
library("VennDiagram")
suppressPackageStartupMessages(library("optparse"))
library("stringr")

option_list <- list(
  make_option(c("--ciri2"), default=NULL, help="File containing results from CIRI2
    identification" ),
  make_option(c("--circexplorer2"), default=NULL, help="File containing results from
    CircExplorer2 identification")
)

parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

if (is.null(opt$ciri2) || is.null(opt$circexplorer2)){
  print_help(parser)
  stop("Results from ciri2 and circexplorer2 identification must be supplied", call.=FALSE)
}else{
  cat(paste("Input files: ", opt$ciri2, " and", opt$circexplorer2))
}

# 1. Load data
ciri_results <- read.csv(opt$ciri2, sep='\t')
circexp_results <- read.csv(opt$circexplorer2, sep='\t', header=FALSE)

# 2. Add circRNA_ID variable in circexplorer2 dataframe
colnames(circexp_results) <- c("chrom", "start", "end", "name", "score", "strand",
  "thickStart", "thickEnd", "itemRgb", "exonCount", "exonSizes", "exonOffsets",
  "readNumber", "circType","geneName","isoformName", "index", "flankIntron")
circRNA_ID <- paste0(circexp_results$chr, ":", circexp_results$start,"|", circexp_results$end)
circexp_results <- cbind(circRNA_ID, circexp_results)

# 3. Filtering ciri2 results by common circRNAs with circexplorer2
circexp_IDs <- circexp_results$circRNA_ID
rownames(ciri_results) <- ciri_results$circRNA_ID
rownames(circexp_results) <- circexp_results$circRNA_ID

coincidences = data.frame()
for (ID in circexp_IDs) {
  if(ID %in% rownames(ciri_results)) {
    if (circexp_results[ID, "strand"] == ciri_results[ID, "strand"]){
      coincidences <- rbind(coincidences, ciri_results[ID, ])
    }
  }
}

# 4. New circRNA matrix
filename= unlist(str_split(opt$ciri2, "libs/identification/ciri2/"))
filename <- filename[2]
filename = unlist(str_split(filename,"_results"))
filename = filename[1]

path_matrix = paste0("libs/identification/", filename, "_coincident_circRNAs.txt")
path_venn = paste0("libs/plots/venn_diagrams/", filename, ".png")

write.table(coincidences, file = path_matrix,
  sep = "\t", row.names = FALSE)

# 5. Grafical representation: Venn Diagramm
venn.diagram(
  x = list(ciri_results$circRNA_ID, circexp_results$circRNA_ID),
  category.names = c("CIRI2", "CircExplorer2"),
  filename = path_venn,
  output = TRUE,

  # Output features
  imagetype="png" ,
  compression = "lzw",

  # Circles
  fill = c("lightblue", "palegreen"),

  # Numbers
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  ext.text = FALSE,

  # Set names
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
)

cat("Completed!!")
paste0("Output files: ", path_matrix, " && ", path_venn)
