#!/bin/R
################################################################################
# Differential expression analysis.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 25-03-2021
# Last modification : 15-06-2021
################################################################################

suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("argparse"))

source("utils/utils.R")

parser <- ArgumentParser(description = 'Differential Expression Analysis')
parser$add_argument("--design", action = "store", type = "character", default = NULL,
              help = "experimental design", metavar = "design")
parser$add_argument("--metadata", type = "character", default = NULL,
            help = "File containing information about samples:Sample names, Group and Sex",
            metavar = "FILE1")
parser$add_argument("--sep", default = ",", action = 'store', help = "Separator")
parser$add_argument("--counts", default = NULL, help = "Numeric matrix of counts",
            metavar="FILE2")
parser$add_argument("--pval", default = 0.05, help = "pvalue", metavar = "pval" )
parser$add_argument("--fc", default = 1.5, help = "foldchange", metavar = "FC" )
parser$add_argument("--outdir", default = NULL, help = "output file")

opt    <- parser$parse_args()

if(is.null(opt$design)|| is.null(opt$metadata) || is.null(opt$counts)){
  parser$print_help()
  stop("Options --design/--metadata/--circ_counts must be supplied\n", call.=FALSE)
}else{
  cat("The supplied arguments are the followings:\n")
  cat(paste(" Experimental design      = ", opt$design, "\n",
            "Library file              = ", opt$metadata, "\n",
            "Count matrix       = ", opt$counts, "\n",
          )
        )
}

# 1. Load data
print("Loading data...")
sep            <- check_sep(opt$sep)
metadata       <- check_metadata(opt$metadata, sep)
counts         <- as.matrix(read.csv(opt$counts, row.names = 1))
design     <- as.formula(opt$design)
print("Done.")

#~~~~~~~~~~~~~~~~~~~~~~~~~DIFFERENTIAL EXPRESSION ANALYSIS~~~~~~~~~~~~~~~~~~~~~~
# 1. DE analysis

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = metadata,
                              design    = design)
dds <- DESeq(dds)
resultsNames(dds)

## Results
# res1 <- results(dds, alpha = 0.05, name = "group_Multiple.sclerosis_vs_Healthy.controls")
# res2 <- results(dds, alpha = 0.05, name = "sex_male_vs_female")
# summary(res1)
# summary(res2)
#
# res1_ordered <- res1[order(res1$padj), ]
# res2_ordered <- res2[order(res2$padj), ]
#
# ## Log FC shrinkage for visualization and ranking
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # The shrunken log fold changes are useful for ranking and visualization,
# # without the need for arbitrary filters on low count genes
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BiocManager::install("apeglm")
# suppressPackageStartupMessages("apeglm")
# res1LFC <- lfcShrink(dds, coef = "group_Multiple.sclerosis_vs_Healthy.controls", type="apeglm")
# res1LFC <- res1LFC[order(res1LFC$padj,)]
# res2LFC <- lfcShrink(dds, coef = "sex_male_vs_female", type="apeglm")
# res2LFC <- res2LFC[order(res2LFC$padj,)]
#
#
# # 2. Exploring and exporting results
# ## MA-plot
# ylim = c(-3 , 3) ; xlim = c (1, 1e5)
# plotMA(res1, xlim= xlim, ylim = ylim, main= "normal")
# plotMA(res1LFC, xlim= xlim, ylim = ylim, main= "apeglm")
# plotMA(res2, xlim= xlim, ylim = ylim, main= "normal")
# plotMA(res2LFC, xlim= xlim, ylim = ylim, main= "apeglm")
#
# ## Plot counts
# plotCounts(dds, gene=which.max(res1$padj), intgroup="group")
# plotCounts(dds, gene=which.max(res2$padj), intgroup="sex")
#
# # Reproducibility
# Sys.time() # date generated
# proc.time() # time spent making
# options(with = 120); sessioninfo::session_info() # R and packages info
