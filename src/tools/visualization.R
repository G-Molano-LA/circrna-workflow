#!/bin/R

###############################################################################
# R script for visualization data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 09-04-2021
# Last modification : 07-05-2021
###############################################################################

# Dependencies
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("factoextra"))

# Functions
check_metadata <- function(metadata){
  if('Group' %in% colnames(metadata) & 'Sex' %in% colnames(metadata) ){
    metadata       <- metadata %>% rename(group = Group)
    metadata       <- metadata %>% rename(sex = Sex)

    metadata$group <- as.factor(metadata$group)
    metadata$sex   <- as.factor(metadata$sex)
    metadata       <- cbind(metadata['group'], metadata['sex'])
    return(metadata)
  }else if('Group' %in% colnames(metadata)){
    metadata            <- metadata %>% rename(group = Group)
    metadata       <- as.factor(metadata['group'])
    return(metadata)
  }else{
    stop(paste0("ERROR: 'Group' column must be supplied in metadata"))
  }
}

check_norm <- function(norm, metadata, circ_info, design){
  design <- as.formula(design)

  if(norm == "False" && circ_info != 'None'){
    circ_info             <- read.csv(opt$circ_info, row.names = 1)

    rownames(circ_counts) <- rownames(circ_info)
    colnames(circ_counts) <- rownames(metadata)
    DESeq_count_data      <- DESeqDataSetFromMatrix(countData = circ_counts,
                                                    colData   = metadata,
                                                    design    = design)
    DESeq_count_data <- estimateSizeFactors(DESeq_count_data)
    circ_counts      <- as.matrix(counts(DESeq_count_data, normalized = TRUE))
    return(circ_counts)
  }else if(norm == "False" && circ_info == 'None'){
    circ_counts      <- as.matrix(read.csv(opt$data, row.names = 1))
    DESeq_count_data <- DESeqDataSetFromMatrix(countData = circ_counts,
                                                colData  = metadata,
                                                design   = design)
    DESeq_count_data <- estimateSizeFactors(DESeq_count_data)
    circ_counts      <- as.matrix(counts(DESeq_count_data, normalized = TRUE))
    return(circ_counts)
  }else if(norm == "True"){
    return(circ_counts)
  }
}

################################################################################
# 1. Input
################################################################################
parser <- ArgumentParser(description = 'Visualization')
parser$add_argument( "--data", default = NULL, help = "Normalized Count Matrix")
parser$add_argument("--norm", default = "no", help = "Boolean variable")
parser$add_argument("--design", default = "~group", help = "Experimental design")
parser$add_argument("--lib", default = NULL, help =
  "File containing library information about samples: sample names, total reads, mapped reads, circular reads, Group and Sex")
parser$add_argument("--circ_info", default = NULL, help =
  "Circular information file, containning circRNA annotation information.")
parser$add_argument( "--output",type = "character", default = "svg")
parser$add_argument("--outdir", type = "character", default = NULL)

opt <- parser$parse_args()

# Load data
metadata    <- read.csv(opt$lib, row.names = 1)
metadata    <- check_metadata(metadata)

circ_counts <- as.matrix(read.csv(opt$data))
circ_counts <- check_norm(opt$norm, metadata, opt$circ_info, opt$design)

################################################################################
# 2. Visualization
################################################################################
circ_counts2 <- melt(circ_counts)
colnames(circ_counts2) <- c("Circular_RNAs", "Samples", "value")
#~~~~~~~~~~~~~~~~~~~~~BOXPLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot
boxplot <-
ggplot(circ_counts2, aes(x = factor(Samples), y = value, fill = factor(Samples))) +
  geom_boxplot() +
  geom_jitter(color="grey", size=0.4, alpha=0.9) +
  labs(fill = "Samples") +
  xlab("Samples") +
  ylab("Normalized counts")

#~~~~~~~~~~~~~~~~~~~~~VIOLINPLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
violinplot <-
ggplot(circ_counts2, aes(x=factor(Samples), y = value, fill = factor(Samples))) +
  geom_violin() +
  labs(fill = "Samples") +
  xlab("Samples") +
  ylab("Normalized counts")

#~~~~~~~~~~~~~~~~~~~~~HISTOGRAM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
histogram <-
ggplot(circ_counts2, aes(x=value)) +
  geom_histogram() +
  ylab("Frequency counts")
  xlab("Number of reads")

#~~~~~~~~~~~~~~~~~~~~~DENDROGRAM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Data treatment
circ_counts3 <- t(circ_counts) # row=samples, columns=circRNAs
dist1        <- dist(circ_counts3, method = "euclidean") # Euclidean distance
dendro       <- hclust(dist1, method = "average") #  Hierarchical Clustering with hclust

# Save Plot
dendro      <- as.dendrogram(dendro)
upper_ylim  <- round(max(dist1), digits = -2) # round to hundredths

svg(filename = paste0(opt$outdir, "/dendrogram.svg"))
plot(dendro, ylab = "Euclidean Distance",xlab = "Samples", main = "Cluster Dendrogram",
  ylim = c(0,upper_ylim))
dev.off()

#~~~~~~~~~~~~~~~~~~~~~PRINCIPAL COMPONENTS ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# more info : http://www.sthda.com/english/wiki/wiki.php?id_contents=7851
# Data treatment
# pca <- prcomp(circ_counts3, center = TRUE, scale = TRUE)
#
# # General plots
# scree_plot         <- fviz_screeplot(pca, ncp = 10)
# graph_of_variables <- fviz_pca_var(pca, col.var = "contrib") +
#   theme_minimal()
# # Individual plots
# dim1_dim2 <- fviz_pca_ind(pca, col.ind="cos2") +
#   scale_color_gradient2(low="blue", mid="white",
#                     high="red", midpoint=0.50)
# dim1_dim3 <- fviz_pca_ind(pca, axes = c(1,3), col.ind="cos2") +
#   scale_color_gradient2(low="blue", mid="white",
#                     high="red", midpoint=0.50)
# dim2_dim3 <- fviz_pca_ind(pca, axes = c(2,3), col.ind="cos2") +
#   scale_color_gradient2(low="blue", mid="white",
#                     high="red", midpoint=0.50)

################################################################################
# 3. Save ggplots
################################################################################

plots <- list(boxplot = boxplot, violin_plot = violinplot,
              histogram = histogram)#, scree_plot = scree_plot,
              #graph_of_variables = graph_of_variables) dim1_dim2 = dim1_dim2,
              #dim1_dim3 = dim1_dim3, dim2_dim3 = dim2_dim3)

for (i in 1:length(plots)){
  switch(opt$output,
    svg = ggsave(filename = paste0(opt$outdir, "/", names(plots[i]), ".svg"),
                 plot = plots[[i]], device = "svg"),
    pdf = ggsave(filename = paste0(opt$outdir, "/", names(plots[i]), ".pdf"),
                 plot = plots[[i]], device = "pdf")
  )
}

switch(opt$output,
       svg = ggsave(filename = paste0(opt$outdir, "/", names(plots[i]), ".svg"),
                    plot = plots[[i]], device = "svg"),
       pdf = ggsave(filename = paste0(opt$outdir, "/", names(plots[i]), ".pdf"),
                    plot = plots[[i]], device = "pdf")
)
