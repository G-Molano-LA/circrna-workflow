#!bin/bash/R

###############################################################################
## Differential expression analysis
###############################################################################
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("statmod"))

# 1. Data
metadata<-read.csv("~/circrna-workflow/docs/metadata_samples.csv",
                   stringsAsFactors=TRUE)
          # data frame containing information for each sample

# shell: wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159225/suppl/GSE159225%5FcircRNA%5FReadcount%5FAllsamples%2Ecsv%2Egz
#         -O docs/circrna_counts.csv
#        wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159225/suppl/GSE159225%5FlinearRNA%5FReadcount%5FAllsamples%2Ecsv%2Egz
#         -O docs/linear_counts.csv
linear_counts <- read.csv("~/circrna-workflow/docs/linear_counts.csv", sep=";")
circrna_counts <- read.csv("~/circrna-workflow/docs/circrna_counts.csv", sep=";")

# Filtering the columns by selected samples
keep <- as.vector(metadata$patient)
print("The selected samples are:"); print(keep)


# 2. The DGEList data class
## Linears
linear_info <- linear_counts[c(1:7)] # data frame containing annotation information
                                    #+ for each gene
linear_counts <- as.matrix(linear_counts[keep])   # numeric matrix of read counts
circrna_info <- circrna_counts[c(1:7)]
circrna_counts <- as.matrix(circrna_counts[keep])

write.csv(linear_info, file = "~/circrna-workflow/libs/DE_analysis/linear_info.csv", quote=FALSE, row.names = FALSE)
write.csv(linear_counts, file = "~/circrna-workflow/libs/DE_analysis/linear_counts_sub.csv", quote=FALSE, row.names = FALSE)
write.csv(circrna_info, file = "~/circrna-workflow/libs/DE_analysis/circular_info.csv", quote=FALSE, row.names = FALSE)
write.csv(circrna_counts, file = "~/circrna-workflow/libs/DE_analysis/circular_counts_sub.csv", quote=FALSE, row.names = FALSE)
