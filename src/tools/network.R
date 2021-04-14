

source = "https://github.com/czllab/NetMiner/blob/master/manual.txt"

# Requirements
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")

bioc_pkg <- c("impute", "preprocessCore", "GO.db")
BiocManager::install(bioc_pkg)

pkg <- c("WGCNA", "corpcor", "bc3net", "igraph", "GeneNet")
install.packages(pkg, dependencies = TRUE)
# Rscript --vanilla ensemble_method_for_construction_consensus_network.R \
#         filePath nThreads percThreshold S1N S2N

# Normalized counts

## 1. Normalized FPKM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The following function returns fragment counts normalized per kilobase of
# feature length per mil-lion mapped fragments (by default using a robust estimate
#  of the library size, as inestimateSizeFactors).
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
counts <- as.matrix(read.csv("libs/DE_analysis/circular_counts_sub.csv"))
linear_counts <- as.matrix(read.csv("libs/DE_analysis/linear_counts_sub.csv"))
circ_info <- read.csv("libs/DE_analysis/circular_info.csv")
linear_info <- read.csv("libs/DE_analysis/linear_info.csv")
basepairs<- circ_info$end-circ_info$start # length of circRNAs
row_names<- circ_info$id
rownames(counts)<- row_names
metadata <- read.csv("libs/DE_analysis/metadata_samples.csv")


suppressPackageStartupMessages(library("DESeq2"))
DESeq_count_data <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ group + sex)
# RAW COUNTS
raw_counts <- counts(DESeq_count_data)
write.table(raw_counts, file = ("libs/network/raw_counts.txt"),
  row.names = TRUE, col.names = NA, quote = FALSE)

# FPKM COUNTS
mcols(DESeq_count_data)$basepairs <- basepairs
fpkm_counts <- fpkm(DESeq_count_data)
write.table(fpkm_counts,file = ("libs/network/fpkm_counts.txt"),
  row.names = TRUE,quote = FALSE)

## 2. Median Normalization
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This function estimates the size factors using the "median ratio method"
# described by Equation 5 inAnders and Huber (2010).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#  MEDIAN COUNTS
DESeq_count_data <- estimateSizeFactors(DESeq_count_data)
cat("The normalization factors for the different sample runs are",
    sizeFactors(DESeq_count_data))
median_counts <- counts(DESeq_count_data,normalized = TRUE)
write.table(median_counts,file = ("libs/network/median_counts.txt"),
  row.names = TRUE,quote = FALSE)

## 2. Variance Stabilized Data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# blind=TRUE should  be  used  for  comparing  samples  in  a  manner  unbiased
#  by  prior  infor-mation  on  samples,  for  example  to  perform  sample  QA
#  (quality  assurance).blind=FALSE should be used for transforming data for
#  downstream analysis,where the full use of the design information should be
#  mad
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# VST COUNTS
dispersions_blind <- estimateDispersions(DESeq_count_data) # estimateDispersions in DESeq pkg
VST_counts <- getVarianceStabilizedData(dispersions_blind)
write.table(VST_counts,file=("libs/network/VST_counts.txt"), row.names = TRUE,
  col.names = NA,quote=FALSE)

## 3. Trimmed Mean Normalization
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The TMM method implements the trimmed mean of M-values method proposed by
# Robinson andOshlack (2010). By default, the M-values are weighted according to
# inverse variances, as computed by the delta method for logarithms of binomial
# random variables. IfrefColumn is unspecified, then the column whose count-per-million
# upper quartile is closest to the mean upper quartile is set as thereference library.
# The TMMwsp method stands for "TMM with singleton pairing".  This is a variant of TMM
# that is intended to perform better for data with a high proportion of zeros.
# In the TMM method, genes that have zero count in either library are ignored when
# comparing pairs of libraries.  In the TMMwsp method, the positive counts from such
# genes are reused to increase the number of features by which the libraries are compared.
# The singleton positive counts are paired up between the libraries in decreasing order of
# size and then a slightly modified TMM method is applied to the re-ordered libraries.
# IfrefColumnis unspecified, then the column with largest sum of square-root counts
# isused as the reference library.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# TMM COUNTS
suppressPackageStartupMessages(library("edgeR"))

linear_DGE <- DGEList(counts   = linear_counts,
                      samples  = metadata,
                      genes    = linear_info)
linear_keep <- filterByExpr(linear_DGE)
linear_DGE <- linear_DGE[linear_keep, keep.lib.sizes=FALSE]
linear_DGE_TMM <- calcNormFactors(linear_DGE)
linear_DGE_UQ <- calcNormFactors(linear_DGE, method = "upperquartile", p = 0.9)

TMM_DGE <- DGEList(counts = counts,
                       samples = metadata,
                       genes= circ_info,
                       lib.size = linear_DGE$samples[, "lib.size"],
                       norm.factors = linear_DGE_TMM$samples[, "norm.factors"])

for(colIndex in seq(1,ncol(TMM_DGE$counts),by=1)){
  TMM_DGE$counts[,colIndex]=(TMM_DGE$counts[,colIndex])/(TMM_DGE$samples$norm.factors[colIndex])
}
write.table(TMM_DGE$counts, file = ("libs/network/TMM_counts.txt"),
  row.names = TRUE, col.names = NA, quote = FALSE)

## 4. Upper Quartile Normalization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (The upperquartile method is the upper-quartile normalization method of Bullard et al (2010),
#+  inwhich the scale factors are calculated from the 75% quantile of the counts
#+ for each library, afterremoving genes that are zero in all libraries. The idea
#+ is generalized here to allow normalization byany quantile of the count distributions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# UQ COUNTS
UQ_DGE <- DGEList(counts = counts,
                       samples = metadata,
                       genes= circ_info,
                       lib.size = linear_DGE$samples[, "lib.size"],
                       norm.factors = linear_DGE_UQ$samples[, "norm.factors"])

for(colIndex in seq(1,ncol(UQ_DGE$counts),by=1)){
  UQ_DGE$counts[,colIndex]=(UQ_DGE$counts[,colIndex])/(UQ_DGE$samples$norm.factors[colIndex])
}

write.table(UQ_DGE$counts, file = ("libs/network/UQ_counts.txt"),
  row.names = TRUE, col.names = NA, quote = FALSE)
