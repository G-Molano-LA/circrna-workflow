#!/bin/R

###############################################################################
## R script for selecting common identified circRNAS by CIRI2 & CircExplorer2
## Author: G. Molano, LA (gonmola@hotmail.es)
###############################################################################
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tibble"))


# 1. Load data
ciri_results <- read.csv("test/SRR12794713_ann_results", sep='\t')
circexp_results <- read.csv("test/SRR12794713_circularRNA_known.txt",
 sep='\t', header=FALSE)

# Reformat 1-based CIRI2 start coordinate to 0-based
ciri_results$circRNA_start <- ciri_results$circRNA_start - 1
ciri_results$circRNA_ID <- paste0(ciri_results$chr, ":", ciri_results$circRNA_start,
  "-", ciri_results$circRNA_end)

# 2. Add circRNA_ID variable in circexplorer2 dataframe
colnames(circexp_results) <- c("chrom", "start", "end", "name", "score", "strand",
  "thickStart", "thickEnd", "itemRgb", "exonCount", "exonSizes", "exonOffsets",
  "readNumber", "circType","geneName","isoformName", "index", "flankIntron")

circexp_results$circRNA_ID <- paste0(circexp_results$chr, ":", circexp_results$start,"-",
  circexp_results$end) # CircExplorer start coordinate is 0-based

# 3. Filtering ciri2 results by common circRNAs with circexplorer2
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
