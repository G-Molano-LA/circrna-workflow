#!/bin/R

################################################################################
# Rscript to normalize raw counts with 5 different methods
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Last modification : 30-06-2021
################################################################################

suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))

source("src/utils/utils.R")

# SCRIPT ARGUMENTS
parser <- ArgumentParser(description = "Normalizing raw counts")
parser$add_argument("--circ_counts", action = 'store', default = NULL, help =
  "circular counts matrix")
parser$add_argument("--metadata", action = 'store', default = NULL, help =
    "Metadata information")
parser$add_argument("--sep", default = ",", action = 'store', help = "Separator")
parser$add_argument("--design", action = 'store', default = NULL, help =
    "Experimental design")
parser$add_argument("--circ_info", action = 'store', default = NULL, help =
    "Circular information")
parser$add_argument("--outdir", action = 'store', default = NULL, help = "output directory")

opt <- parser$parse_args()

if(is.null(opt$circ_counts)){
  parser$print_help()
  stop("Options --circ_counts/--metadata must be supplied\n", call. = FALSE)
}


# LOAD DATA
sep         <- check_sep(opt$sep)
metadata    <- check_metadata(opt$metadata, sep)
design      <- as.formula(opt$design)
circ_info   <- read.csv(opt$circ_info)
circ_counts <- check_data_co(opt$circ_counts, circ_info)

DESeq_count_data <- DESeqDataSetFromMatrix(countData = circ_counts,
                                            colData  = metadata,
                                            design = design)
DGE              <- DGEList(counts = circ_counts,
                            samples= metadata)

# 1. Raw Counts
raw_counts <- counts(DESeq_count_data)


# 2. FPKM Counts
# 2.1. Calculate circRNA length. Only exoncircRNAs.
list <- strsplit(circ_info$exon_sizes, ",")
circ_info <- add_column(circ_info, length = NA, .after = "exon_sizes")
for (i in 1:length(list)) {
  circ_info$length[i] <- sum(as.numeric(list[[i]]))
}
mcols(DESeq_count_data)$basepairs <- circ_info$length # length of circRNAs
fpkm_counts <- fpkm(DESeq_count_data)


# 3. Median Counts
DESeq_count_data <- estimateSizeFactors(DESeq_count_data)
cat("The normalization factors for the different sample runs are",
    sizeFactors(DESeq_count_data))

median_counts <- counts(DESeq_count_data,normalized = TRUE)


# 4. VST Counts (Variance Stabilized Data)
dispersions_blind <- estimateDispersions(DESeq_count_data) # estimateDispersions in DESeq pkg
vst_counts        <- getVarianceStabilizedData(dispersions_blind)


# 5. Trimmed Mean Counts (TMM)
DGE_TMM <- calcNormFactors(DGE)

for(colIndex in seq(1,ncol(DGE_TMM$counts),by=1)){
  DGE_TMM$counts[,colIndex]=(DGE_TMM$counts[,colIndex])/(DGE_TMM$samples$norm.factors[colIndex])
}


# 6. Upper Quartile Counts (UQ)
DGE_UQ <- calcNormFactors(DGE, method = "upperquartile", p = 0.9)

for(colIndex in seq(1,ncol(DGE_UQ$counts),by=1)){
  DGE_UQ$counts[,colIndex]=(DGE_UQ$counts[,colIndex])/(DGE_UQ$samples$norm.factors[colIndex])
}

# DOWNLOAD
files <- list(Raw = raw_counts, FPKM = fpkm_counts, Median = median_counts,
              VST = vst_counts, TMM = DGE_TMM$counts, UQ = DGE_UQ$counts)

for(i in 1:length(files)){
  write.table(files[[i]],
              file = paste0(opt$outdir, names(files[i]), "_Count.txt"),
              row.names = TRUE, quote = FALSE, sep = "\t")
}
