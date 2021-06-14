#!/bin/R

# 1. DE analysis
counts <- as.matrix(read.csv("libs/DE_analysis/circular_counts_sub.csv"))
circ_info <- read.csv("libs/DE_analysis/circular_info.csv")
row_names<- circ_info$id
rownames(counts)<- row_names
metadata <- read.csv("libs/DE_analysis/metadata_samples.csv")


suppressPackageStartupMessages(library("DESeq2"))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ group + sex)
dds <- DESeq(dds)
resultsNames(dds)

## Results
res1 <- results(dds, alpha = 0.05, name = "group_Multiple.sclerosis_vs_Healthy.controls")
res2 <- results(dds, alpha = 0.05, name = "sex_male_vs_female")
summary(res1)
summary(res2)

res1_ordered <- res1[order(res1$padj), ]
res2_ordered <- res2[order(res2$padj), ]

## Log FC shrinkage for visualization and ranking
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The shrunken log fold changes are useful for ranking and visualization,
# without the need for arbitrary filters on low count genes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BiocManager::install("apeglm")
suppressPackageStartupMessages("apeglm")
res1LFC <- lfcShrink(dds, coef = "group_Multiple.sclerosis_vs_Healthy.controls", type="apeglm")
res1LFC <- res1LFC[order(res1LFC$padj,)]
res2LFC <- lfcShrink(dds, coef = "sex_male_vs_female", type="apeglm")
res2LFC <- res2LFC[order(res2LFC$padj,)]


# 2. Exploring and exporting results
## MA-plot
ylim = c(-3 , 3) ; xlim = c (1, 1e5)
plotMA(res1, xlim= xlim, ylim = ylim, main= "normal")
plotMA(res1LFC, xlim= xlim, ylim = ylim, main= "apeglm")
plotMA(res2, xlim= xlim, ylim = ylim, main= "normal")
plotMA(res2LFC, xlim= xlim, ylim = ylim, main= "apeglm")

## Plot counts
plotCounts(dds, gene=which.max(res1$padj), intgroup="group")
plotCounts(dds, gene=which.max(res2$padj), intgroup="sex")

# Reproducibility
Sys.time() # date generated
proc.time() # time spent making
options(with = 120); sessioninfo::session_info() # R and packages info
