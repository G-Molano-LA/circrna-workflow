#!/bin/R

###############################################################################
## R script for visualizing results from CIRI2 & CircExplorer2, as well as,
## common identified circRNAs by overlap function.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date : 30-06-2021
################################################################################

suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tibble"))

parser <- ArgumentParser(description = 'circRNA identified list')
parser$add_argument("--ciri2", default = NULL, help = "File containing results from CIRI2 identification.")
parser$add_argument("--circexplorer2", default = NULL, help = "File containing results from CircExplorer2 identification.")
parser$add_argument("--overlap", default = NULL, help = "File containing common identified circRNAs from both tools.")
parser$add_argument("--outdir", action = 'store', default = NULL, help = "output directory")

opt <- parser$parse_args()

if (is.null(opt$ciri2) || is.null(opt$circexplorer2) || is.null(opt$outdir)){
  print_help(parser)
  stop("Results from ciri2 and circexplorer2 identification must be supplied", call.=FALSE)
}else{
  cat(paste("Input files: ", opt$ciri2, " and", opt$circexplorer2))
}

# 1. Load data
ciri_results    <- read.csv(opt$ciri2, sep='\t')
circexp_results <- read.csv(opt$circexplorer2, sep='\t')
overlap_results <- read.csv(opt$overlap, sep = '\t')


# 6. New circRNA matrix
path_venn <- paste0(opt$outdir, "summary_overlap.png")

# 7. Grafical representation: Venn Diagramm
venn <- draw.pairwise.venn(
  area1 = nrow(circexp_results),
  area2 = nrow(ciri_results),
  cross.area = nrow(overlap_results),
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
paste0("Output files: ", path_venn)
