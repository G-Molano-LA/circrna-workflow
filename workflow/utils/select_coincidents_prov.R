#!/bin/R

###############################################################################
## R script for selecting common identified circRNAS by CIRI2 & CircExplorer2
## Author: G. Molano, LA (gonmola@hotmail.es)
###############################################################################
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stringr"))


# 1. Load data
ciri_results <- read.csv("libs/identification/ciri2/SRR12794713_ann_results", sep='\t')
circexp_results <- read.csv("libs/identification/circexplorer2/SRR12794713_circularRNA_known.txt", sep='\t', header=FALSE)

# Reformat ciri2 circRNA_ID variable = to be more flexible
start2 <- substr(ciri_results$circRNA_start, 1, nchar(ciri_results$circRNA_start)-1)
end2 <- substr(ciri_results$circRNA_end, 1, nchar(ciri_results$circRNA_end)-1)
ciri_results$circRNA_ID <- paste0(ciri_results$chr, ":", start2,"|", end2)


# 2. Add circRNA_ID variable in circexplorer2 dataframe
colnames(circexp_results) <- c("chrom", "start", "end", "name", "score", "strand",
  "thickStart", "thickEnd", "itemRgb", "exonCount", "exonSizes", "exonOffsets",
  "readNumber", "circType","geneName","isoformName", "index", "flankIntron")

start1 <- substr(circexp_results$start, 1, nchar(circexp_results$start)-1)
end1 <- substr(circexp_results$end, 1, nchar(circexp_results$end)-1)
circRNA_ID <- paste0(circexp_results$chr, ":", start1,"|", end1)
circexp_results <- cbind(circRNA_ID, circexp_results)

# 3. Filtering ciri2 results by common circRNAs with circexplorer2
coincidences = data.frame()
for (id in circexp_results$circRNA_ID) {
  if(id %in% ciri_results$circRNA_ID) {
    index <- which(ciri_results$circRNA_ID == id)
    coincidences <- rbind(coincidences, ciri_results[index, ])
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
  filename = "risas.png",
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
