#!/bin/R

NOTA = "Falta modificar en funció dels inputs que agafi. Enteoria, agafarè... pos no se jajajajja "
###############################################################################
## R script for annotating circular RNAs with circBase ID
## Author: G. Molano, LA
###############################################################################
source = "https://bioconductor.org/packages/release/workflows/vignettes/liftOver/inst/doc/liftov.html#resource-the-chain-file-for-hg38-to-hg19-transformation"

suppressPackageStartupMessages(library(rtracklayer))

# Load data
circbase <- read.csv("~/Downloads/hsa_hg19_circRNA.txt", sep = "\t", header = TRUE)
list <- read.csv("libs/identification/ciri2/SRR12794713_ann_results", sep="\t",
  header = TRUE)

# Reformat 1-based CIRI2 start coordinate to 0-based: enteoria això ja estarà fet.
list$circRNA_start <- list$circRNA_start - 1
list$circRNA_ID <- paste0(list$chr, ":", list$circRNA_start,
  "-", list$circRNA_end)

# The chain file provided by UCSC for hg38 to hg19 transformation
# wget --timestamping
#         'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz'
#         -O hg38ToHg19.over.chain.gz
chain = import.chain("data/raw_data/hg38ToHg19.over.chain")

# Convert dataframe into GRanges list
list_range <- makeGRangesFromDataFrame(
  df = list,
  start.field="circRNA_start",
  end.field="circRNA_end",
  keep.extra.columns = TRUE
)

# liftOver conversion
seqlevelsStyle(list_range) = "UCSC"
list19 <- liftOver(list_range, chain)
list19 <- unlist(list19)
genome(list19) <- "hg19"

# Convert GRanges list into dataframe
list19 <- as.data.frame(list19)
list19 <- subset(list19, select = - circRNA_ID)
list19$circRNA_ID <- paste0(list19$seqnames, ":", list19$start,
  "-", list19$end)

# Find circBase ID
circbase$circRNA_ID <- paste0(circbase$X..chrom, ":", circbase$start,
  "-", circbase$end)

list19$circBase_ID <- NA
for (id in list19$circRNA_ID) {
  if(id %in% circbase$circRNA_ID) {
    index1 <- which(circbase$circRNA_ID == id)
    index2 <- which(list19$circRNA_ID == id)

    circBase_ID <- circbase$circRNA.ID[index1]
    list19$circBase_ID[index2] <- circBase_ID
  }
}

# Save data
write.table(list19, file = "annotated_circRNAs.txt",sep = "\t", row.names = FALSE)
