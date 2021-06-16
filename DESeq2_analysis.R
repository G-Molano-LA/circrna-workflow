suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("apeglm"))

source("~/circrna-workflow/workflow/utils/utils.R")

sep            <- check_sep("comma")
metadata       <- check_metadata("~/circrna-workflow/libs/DE_analysis/metadata.csv", sep)
counts         <- as.matrix(read.csv("~/circrna-workflow/results/circ_counts.txt", sep = "\t", row.names = 1))
design         <- as.formula("~group+sex")

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = metadata,
                              design    = design)
dds <- DESeq(dds)
resultsNames(dds)

res1 <- results(dds, alpha = 0.05, name = "group_Multiple.sclerosis_vs_Healthy.controls")
res2 <- results(dds, alpha = 0.05, name = "sex_male_vs_female")

summary(res1)
summary(res2)

res1_ordered <- res1[order(res1$padj), ]
res2_ordered <- res2[order(res2$padj), ]

DE_group <- sum(res1_ordered$padj < 0.05, na.rm = TRUE)
DE_sex   <- sum(res2_ordered$padj < 0.05, na.rm = TRUE)

print(paste("Group Multiple Sclerosis vs Healthy controls: ", DE_group, " differentially expressed circRNAs." ))
print(paste("Sex male vs female: ", DE_sex, " differentially expressed circRNAs." ))

ylim = c(-3 , 3) ; xlim = c (1, 1e3)
plotMA(res1, xlim = xlim, ylim = ylim, main = "normal", colSig = "red", alpha = 0.05)
plotMA(res2, xlim = xlim, ylim = ylim, main = "normal", colSig = "red", alpha = 0.05)
