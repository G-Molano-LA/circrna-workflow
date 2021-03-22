#!bin/bash/R

###############################################################################
## Differential expression analysis
###############################################################################
library(edgeR)

# 1. Data
metadata<- read.csv("~/Desktop/TFG/data/esclerosis/processed_data/metadata_sampled.csv")
          # data frame containing information for each sample
counts <- read.csv("~/Desktop/TFG/data/esclerosis/processed_data/GSE159225_circRNA_Readcount_Allsamples.csv", sep=";")

# Filtering the columns by selected samples
keep <- metadata$patient
print("The selected samples are:"); print(keep)


# 2. The DGEList data class
circrna_counts <- counts[keep]
as.matrix(circrna_counts)      # numeric matrix of read counts
circrna_info <- counts[c(1:7)] # data frame containing annotation information 
                               #+ for each circrna

circ_DGE <- DGEList(counts=circrna_counts, samples=metadata, 
             group=metadata$disease_state, genes=circrna_info, remove.zeros = T)
