#!bin/bash/R

###############################################################################
## Differential expression analysis
###############################################################################
library("edgeR")
library("statmod")

# 1. Data
metadata<-read.csv("~/alejandra/circrna-workflow/docs/metadata_samples.csv",
                   stringsAsFactors=TRUE)
          # data frame containing information for each sample

# shell: wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159225/suppl/GSE159225%5FcircRNA%5FReadcount%5FAllsamples%2Ecsv%2Egz 
#         -O docs/circrna_counts.csv
#        wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159225/suppl/GSE159225%5FlinearRNA%5FReadcount%5FAllsamples%2Ecsv%2Egz 
#         -O docs/linear_counts.csv
linear_counts <- read.csv("~/alejandra/circrna-workflow/docs/linear_counts.csv", sep=";")
circrna_counts <- read.csv("~/alejandra/circrna-workflow/docs/circrna_counts.csv", sep=";")

# Filtering the columns by selected samples
keep <- metadata$patient
keep <- as.vector(keep)
print("The selected samples are:"); print(keep)


# 2. The DGEList data class
## Linears
linear_info <- linear_counts[c(1:7)] # data frame containing annotation information 
                                    #+ for each gene
linear_counts <- linear_counts[keep]
linear_counts <- as.matrix(linear_counts)      # numeric matrix of read counts


linear_DGE <- DGEList(counts   = linear_counts, 
                      samples  = metadata, 
                      group    = metadata$disease_state, 
                      genes    = linear_info, remove.zeros = T)
linear_idx <- filterByExpr(linear_DGE)
linear_DGE <- linear_DGE[linear_idx, keep.lib.sizes=FALSE]
linear_DGE <- calcNormFactors(linear_DGE)
design <- model.matrix (~metadata$disease_state + metadata$sex)

## Circulars
circrna_info <- circrna_counts[c(1:7)]
circrna_counts <- circrna_counts[keep] 
circrna_DGE <- DGEList(counts = circrna_counts,
                       samples = metadata,
                       group = metadata$disease_state,
                       genes= circrna_info,
                       lib.size = linear_DGE$samples[, "lib.size"],
                       norm.factors = linear_DGE$samples[, "norm.factors"])

circrna_DGE <- estimateDisp(circrna_DGE, design, robust = TRUE)
circrna_fit <- glmFit(circrna_DGE, design)
circrna_lrt <- glmLRT(circrna_fit)

circrna_df <- circrna_lrt$table
circrna_order <- order(circrna_lrt$table$PValue)
circrna_df$DE <- as.vector(decideTestsDGE(circrna_lrt))
circrna_df <- circrna_df[circrna_order, ]
circrna_df$FDR <- p.adjust(circrna_df$PValue, method="fdr")
write.csv(circrna_df, file = "~/alejandra/circrna-workflow/docs/circrna_DE.csv", quote=FALSE)
