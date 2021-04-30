#!/bin/R

###############################################################################
## R script for selecting common identified circRNAS by CIRI2 & CircExplorer2
## Author: G. Molano, LA (gonmola@hotmail.es)
###############################################################################
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tibble"))

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
ciri_results    <- read.csv(opt$ciri2, sep='\t')
circexp_results <- read.csv(opt$circexplorer2, sep='\t', header=FALSE)

path <- unlist(str_split(opt$ciri2, "/identification/")) # split the beggining
outdir <- path[1]
path <- path[2]
path <- unlist(str_split(path,"_ciri2.txt")) # split the end
sample <- path[1]

path_matrix <- paste0(outdir, "/identification/overlap/", sample, "_common.txt")
path_venn   <- paste0(outdir, "/identification/overlap", sample, "_common.png")

# 2. Reformat 1-based CIRI2 start coordinate to 0-based
ciri_results$circRNA_start <- ciri_results$circRNA_start - 1
ciri_results$circRNA_ID    <- paste0(ciri_results$chr, ":", ciri_results$circRNA_start,
  "-", ciri_results$circRNA_end)

# 3. Add circRNA_ID variable in circexplorer2 dataframe
colnames(circexp_results) <- c("chrom", "start", "end", "name", "score", "strand",
  "thickStart", "thickEnd", "itemRgb", "exonCount", "exonSizes", "exonOffsets",
  "readNumber", "circType","geneName","isoformName", "index", "flankIntron")

circexp_results$circRNA_ID <- paste0(circexp_results$chr, ":", circexp_results$start,"-",
  circexp_results$end) # CircExplorer start coordinate is 0-based

# 4. Filtering ciri2 results by common circRNAs with circexplorer2
overlap = data.frame()
ciri_results <- add_column(ciri_results, gene_id = NA, .after = "circRNA_type") # adding genename from circexplorer2 results

for (id in circexp_results$circRNA_ID) {
  if(id %in% ciri_results$circRNA_ID) {
    index1 <- which(circexp_results$circRNA_ID == id)
    index2 <- which(ciri_results$circRNA_ID == id)

    ciri_results$gene_id[index2] <- circexp_results$geneName[index1]
    overlap <- rbind(overlap, ciri_results[index2, ])
  }
}
overlap <- overlap[!duplicated(overlap), ]


# 5. New circRNA matrix

write.table(overlap, file = path_matrix, quote = FALSE, sep = "\t", row.names = FALSE)

# 6. Grafical representation: Venn Diagramm
venn <- draw.pairwise.venn(
  area1 = nrow(circexp_results),
  area2 = nrow(ciri_results),
  cross.area = nrow(overlap),
  category = c("CircExplorer2", "CIRI2"),

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
