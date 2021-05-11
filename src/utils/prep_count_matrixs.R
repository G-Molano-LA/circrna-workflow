#!bin/R

###############################################################################
# R script for preparing inputs for DE analysis. Initial data from ncbi repository
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 
# Last modification : 10-05-2021
################################################################################
suppressPackageStartupMessages(library("statmod"))
suppressPackageStartupMessages(library("dplyr"))

setwd("~/circrna-workflow")

# 1. METADATA file
metadata<-read.csv("~/circrna-workflow/docs/metadata_all.txt",
                   stringsAsFactors=TRUE) # Download from SRA repository
metadata_filtered <- subset(metadata,  select=c(Run, disease_state, ms_type, sex))
metadata_filtered <- metadata_filtered %>% rename(group = disease_state) # must be a colunm name group for edgeR analysis
metadata_filtered <- metadata_filtered[order(metadata_filtered$Run),]

# 1.1. Adding patient id
rrms <- character()
rrms <- c(paste0("RRMS0", 1:9), paste0("RRMS", 10:20))

spms <- character()
spms <- c(paste0("SPMS0", 1:9), paste0("SPMS", 10))

hc <- character()
hc <- c(paste0("HC0", 1:9), paste0("HC", 10:20))


patient_id <- c(rrms, spms, hc)
metadata_filtered$patient <- patient_id

# 1.2. Chosing a subset of samples
# set.seed(100)
# metadata_sample <- metadata_filtered %>% group_by(group, sex) %>% slice_sample(n=2)
# metadata_sample <- metadata_sample[order(metadata_sample$Run), ]
# print("The selected samples are the following patients: ")
# print(metadata_sample$patient)

# 1.3. Writing metadata associated, which contains information for each sample
# write.csv(metadata_sample, file="~/circrna-workflow/libs/DE_analysis/metadata_samples.csv", row.names=F)


# 2. COUNTS file
# system("wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159225/suppl/GSE159225%5FcircRNA%5FReadcount%5FAllsamples%2Ecsv%2Egz \
#   -O docs/circrna_counts.csv", intern = TRUE)
# system("wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159225/suppl/GSE159225%5FlinearRNA%5FReadcount%5FAllsamples%2Ecsv%2Egz \
#   -O docs/linear_counts.csv", intern = TRUE)

circrna_counts <- read.csv("~/circrna-workflow/docs/circrna_counts.csv", sep=";")
linear_counts <- read.csv("~/circrna-workflow/docs/linear_counts.csv", sep=";")


# 2.1. Filtering count matrixs to obtain only counts from selected samples
# keep <- as.vector(metadata_sample$patient)
# print("The selected samples are:"); print(keep)

linear_info <- linear_counts[c(1:7)] # gene annotation information
circrna_info <- circrna_counts[c(1:7)]

# linear_counts <- as.matrix(linear_counts[keep])  # filtered numeric matrix of read counts
# circrna_counts <- as.matrix(circrna_counts[keep])

linear_counts  <- as.matrix(linear_counts[c(8:ncol(linear_counts))])  # filtered numeric matrix of read counts
circrna_counts <- as.matrix(circrna_counts[c(8:ncol(circrna_counts))])

write.csv(linear_info, file = "~/circrna-workflow/libs/DE_analysis/linear_info.csv", quote=FALSE, row.names = FALSE)
write.csv(linear_counts, file = "~/circrna-workflow/libs/DE_analysis/linear_counts_sub.csv", quote=FALSE, row.names = FALSE)
write.csv(circrna_info, file = "~/circrna-workflow/libs/DE_analysis/circular_info.csv", quote=FALSE, row.names = FALSE)
write.csv(circrna_counts, file = "~/circrna-workflow/libs/DE_analysis/circular_counts_sub.csv", quote=FALSE, row.names = FALSE)
write.csv(metadata_filtered, file = "~/circrna-workflow/libs/DE_analysis/metadata.csv", quote=FALSE, row.names = FALSE)
