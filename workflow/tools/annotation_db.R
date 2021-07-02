#!/bin/R

###############################################################################
## R script for annotating circular RNAs with circBase ID
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Last modification : 30-06-2021
###############################################################################

source = "https://bioconductor.org/packages/release/workflows/vignettes/liftOver/inst/doc/liftov.html#resource-the-chain-file-for-hg38-to-hg19-transformation"

suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("dplyr"))

################################################################################
# 1. Input
################################################################################
parser <- ArgumentParser(description = 'circBase_ID Annotation')
parser$add_argument( "--list", default = NULL, help = "File cointainning circRNA list with genomic coordinates (chrom, start, end)")
parser$add_argument( "--db", default = NULL, help = "File cointainning circRNA data from circBase database.")
parser$add_argument( "--chain", default = NULL, help = "File to transform hg38 to hg19 coordinates.")
parser$add_argument("--fields", default = NULL, help = "Extra fields to add.", nargs = '+')
parser$add_argument("--out", type = "character", default = NULL)

opt <- parser$parse_args()

# Load data
list      <- read.csv(opt$list, sep = "\t", header = TRUE)
circbase  <- read.csv(opt$db, sep = "\t", header = TRUE)
chain     <- import.chain(opt$chain) # The chain file provided by UCSC for hg38 to hg19 transformation
fields    <- opt$fields


# CIRI2 start coordinate is 1-based. UCSC annotation is also 1-based
list$circRNA_ID <- paste0(list$X.chrom, ":", list$start, "-", list$end)

# Convert dataframe into GRanges list
list_range <- makeGRangesFromDataFrame(
  df                      = list,
  seqnames.field          = "X.chrom",
  start.field             = "start",
  end.field               = "end",
  strand.field            = "strand",
  starts.in.df.are.0based = FALSE,
  keep.extra.columns      = TRUE
)

# liftOver conversion
seqlevelsStyle(list_range) <- "UCSC"
list19                     <- liftOver(list_range, chain)
list19                     <- unlist(list19)
genome(list19)             <- "hg19"

# Convert GRanges list into dataframe
list19            <- as.data.frame(list19)
list19            <- subset(list19, select = - circRNA_ID)
list19            <- add_column(list19, "circRNA_ID" = NA, .before = "seqnames")
list19$circRNA_ID <- paste0(list19$seqnames, ":", list19$start, "-", list19$end)

# Find circBase ID
circbase$circRNA_ID <- paste0(circbase$X..chrom, ":", circbase$start,
  "-", circbase$end)

list19 <- add_column(list19, "circBase_ID" = NA, .after = "score")
if(fields != "None"){
  list19[fields]     <- NA # creating new columns depending on selected fields

  names <- colnames(list19)[1:9]
  list19 <- list19 %>% select(all_of(names), all_of(fields), everything()) # rearrenging columns
}


for (id in list19$circRNA_ID) {
  if(id %in% circbase$circRNA_ID) {
    index1 <- which(circbase$circRNA_ID == id)
    index2 <- which(list19$circRNA_ID == id)
    if(fields != "None"){
      list19$circBase_ID[index2]<- circbase$circRNA.ID[index1]
      list19[index2, fields]    <- circbase[index1, fields]
    }else{
      list19$circBase_ID[index2]<- circbase$circRNA.ID[index1]
    }
  }
}

# Save data
write.table(list19, file = opt$out ,sep = "\t", row.names = FALSE, quote = FALSE)
