#!/bin/R

###############################################################################
## R script for selecting common identified circRNAS by CIRI2 & CircExplorer2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 16-04-2021
# Last modification : 26-05-2021
################################################################################
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tibble"))

parser <- ArgumentParser(description = 'circRNA identified list')
parser$add_argument("--ciri2", default = NULL, help = "File containing results from CIRI2 identification")
parser$add_argument("--circexplorer2", default = NULL, help = "File containing results from CircExplorer2 identification")
parser$add_argument("--sample", default = NULL, help = "Sample name")
parser$add_argument("--outdir", action = 'store', default = NULL, help = "output directory")

opt <- parser$parse_args()

if (is.null(opt$ciri2) || is.null(opt$circexplorer2) || is.null(opt$sample) || is.null(opt$outdir)){
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

# 3. Reformat 0-based CircExplorer2 start coordinate to 1-based
circexp_results$start <- circexp_results$start + 1
circexp_results$circRNA_ID <- paste0(circexp_results$chr, ":", circexp_results$start,"|",
  circexp_results$end) # CIRI start coordinate is 0-based

# 4. Filtering ciri2 results by common circRNAs with circexplorer2
overlap = data.frame()
ciri_results <- add_column(ciri_results, gene_name = NA, .after = "junction_reads_ID") # adding genename from circexplorer2 results
ciri_results <- add_column(ciri_results, exon_count = NA, .after = "gene_name")
ciri_results <- add_column(ciri_results, exon_sizes = NA, .after = "exon_count")

for (id in circexp_results$circRNA_ID) {
  if(id %in% ciri_results$circRNA_ID) {
    index1 <- which(circexp_results$circRNA_ID == id)
    index2 <- which(ciri_results$circRNA_ID == id)

    ciri_results$gene_name[index2] <- circexp_results$geneName[index1]
    ciri_results$exon_count[index2] <- circexp_results$exonCount[index1]
    ciri_results$exon_sizes[index2] <- circexp_results$exonSizes[index1]
    overlap <- rbind(overlap, ciri_results[index2, ])
  }
}
overlap <- overlap[!duplicated(overlap), ]

# 6. New circRNA matrix
path_matrix = paste0(opt$outdir, opt$sample, "_common.txt")
path_venn = paste0(opt$outdir, opt$sample, "_common.png")

write.table(overlap, file = path_matrix, quote = FALSE, sep = "\t", row.names = FALSE)

# 7. Grafical representation: Venn Diagramm
venn <- draw.pairwise.venn(
  area1 = nrow(circexp_results),
  area2 = nrow(ciri_results),
  cross.area = nrow(overlap),
  category = c(paste0("CircExplorer2\n","(n = ", nrow(circexp_results), ")"),
               paste0("CIRI2\n","(n = ", nrow(ciri_results), ")")),

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
png(filename = path_venn)
grid.draw(venn);
dev.off()


cat("Completed!!")
paste0("Output files: ", path_matrix, " && ", path_venn)
