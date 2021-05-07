#!/bin/R

################################################################################
# Rscript to normalize raw counts with 5 different methods
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 12-04-2021
# Last modification : 06-05-2021
################################################################################

suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("argparse"))

# FUNCTIONS

# !! Hacer una nota al usuario de poner el nombre en de las columnas con la primera
# en mayúscula, siempre y cuando ponga él sus propios archivos
check_lib <- function(lib){
  if(ncol(lib) == 5 ){
    lib      <- lib %>% rename(Group = group)
    metadata <- as.factor(lib['group'])
    return(metadata)
  }elif(ncol(lib) == 6){
    lib            <- lib %>% rename(Group = group)
    lib            <- lib %>% rename(Sex = sex)
    metadata       <- cbind(lib['group'], lib['sex'])
    metadata$group <- as.factor(metadata$group)
    metadata$sex   <- as.factor(metadata$sex)
    return(metadata)
  }else{
    stop(paste0("ERROR: Invalid number of columns in library information file."))
  }
}

# SCRIPT ARGUMENTS
parser <- ArgumentParser(description = "Normalizing raw counts")
parser$add_argument("--circ_counts", action = 'store', default = NULL, help =
  "circular counts matrix")
parser$add_argument("--metadata", action = 'store', default = NULL, help =
    "Metadata information")
parser$add_argument("--circ_info", action = 'store', default = NULL, help =
    "Circular information")
parser$add_argument("--outdir", action = 'store', default = NULL, help = "output directory")

opt <- parser$parse_args()

if(is.null(opt$circ_counts)){
  parser$print_help()
  stop("Options --circ_counts/--metadata must be supplied\n", call. = FALSE)
}


# LOAD DATA
circ_info   <- read.csv(opt$circ_info)

circ_counts           <- as.matrix(read.csv(opt$circ_counts))
rownames(circ_counts) <- circ_info$id

metadata    <- read.csv(opt$metadata, row.names = 1)
metadata    <- check_lib(metadata)

DESeq_count_data <- DESeqDataSetFromMatrix(countData = circ_counts,
                                            colData  = metadata)
DGE              <- DGEList(counts = counts,
                            samples= metadata,
                            genes  = circ_info)

# 1. Raw Counts
raw_counts <- counts(DESeq_count_data)


# 2. FPKM Counts
mcols(DESeq_count_data)$basepairs <- circ_info$end-circ_info$start # length of circRNAs
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
              file = paste0(opt$outdir, names(files[i]), "_Count.txt" )
              row.names = TRUE, quote = FALSE)
}
