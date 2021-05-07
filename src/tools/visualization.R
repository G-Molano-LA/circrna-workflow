#!/bin/R

###############################################################################
# R script for visualization data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 09-04-2021
# Last modification : 06-05-2021
###############################################################################

# Dependencies
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("factoextra"))

check_norm <- function(norm, lib, circ_info){
  if(norm == "no" && circ_info != 'None'){
    circ_info             <- read.csv(opt$circ_info, row.names = 1)
    rownames(circ_counts) <- rownames(circ_info)
    colnames(circ_counts) <- rownames(lib)
    DESeq_count_data      <- DESeqDataSetFromMatrix(countData = circ_counts,
                                                colData  = metadata)
    circ_counts <- counts(DESeq_count_data, normalized = TRUE)
    return(circ_counts)
  }else if(norm == "no" && circ_info == 'None'){
    DESeq_count_data <- DESeqDataSetFromMatrix(countData = circ_counts,
                                                colData  = metadata)
    circ_counts <- counts(DESeq_count_data, normalized = TRUE)
    return(circ_counts)
  }else if(norm == "yes"){
    return(circ_counts)
  }
}

################################################################################
# 1. Input
################################################################################
parser <- ArgumentParser(description = 'Visualization')
parser$add_argument( "--data", default = NULL, help = "Normalized Count Matrix")
parser$add_argument("--norm", default = "no", help = "Boolean variable")
parser$add_argument("--lib", default = NULL, help =
  "File containing library information about samples: sample names, total reads, mapped reads, circular reads, Group and Sex")
parser$add_argument("--circ_info", default = NULL, help =
  "Circular information file, containning circRNA annotation information.")
parser$add_argument( "--output",type = "character", default = "svg")
parser$add_argument("--outdir", type = "character", default = NULL)

opt <- parser$parse_args()

# Load data
circ_counts <- as.matrix(read.csv(opt$data))
metadata    <- read.csv(opt$lib, row.names = 1

circ_counts <- check_norm(opt$norm, metadata, opt$circ_info)

################################################################################
# 2. Visualization
################################################################################
#~~~~~~~~~~~~~~~~~~~~~HEATMAP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Data treatment: Plot top10 DE (by pvalue) circular RNAs
circ_counts1 <- filter(circ_counts, rownames(circ_counts) %in% top20)
circ_counts1 <- as.matrix(circ_counts1)
circ_counts2 <- melt(circ_counts1) # changing format
colnames(circ_counts2) <- c("Circular_RNAs", "Samples", "value")

# Plot
heatmap <-
ggplot(circ_counts2, aes(x=Samples, y=Circular_RNAs, fill=value)) +
  geom_tile(colour = "white") +
  labs(fill = "TMM counts") +
  scale_fill_gradient(low="white", high = "steelblue")

#~~~~~~~~~~~~~~~~~~~~~BOXPLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot
boxplot <-
ggplot(circ_counts2, aes(x = factor(Samples), y = value, fill = factor(Samples))) +
  geom_boxplot() +
  geom_jitter(color="grey", size=0.4, alpha=0.9) +
  labs(fill = "Samples") +
  xlab("Samples") +
  ylab("Normalized TMM counts")

#~~~~~~~~~~~~~~~~~~~~~VIOLINPLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
violinplot <-
ggplot(circ_counts2, aes(x=factor(Samples), y = value, fill = factor(Samples))) +
  geom_violin() +
  labs(fill = "Samples") +
  xlab("Samples") +
  ylab("Normalized TMM counts")

#~~~~~~~~~~~~~~~~~~~~~HISTOGRAM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
histogram <-
ggplot(circ_counts2, aes(x=value)) +
  geom_histogram() +
  xlab("Number of reads")

#~~~~~~~~~~~~~~~~~~~~~DENDROGRAM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Data treatment
circ_counts3 <- t(circ_counts) # row=samples, columns=circRNAs
dist1        <- dist(circ_counts3, method = "euclidean") # Euclidean distance
dendro       <- hclust(dist1, method = "average") #  Hierarchical Clustering with hclust

# Save Plot
dendro      <- as.dendrogram(dendro)
upper_ylim  <- round(max(dist1), digits = -2) # round to hundredths

svg(file = "/libs/plots/dendrograma.svg")
plot(dendro, ylab = "Euclidean Distance",xlab = "Samples", main = "Cluster Dendrogram",
  ylim = c(0,upper_ylim))
dev.off()

#~~~~~~~~~~~~~~~~~~~~~PRINCIPAL COMPONENTS ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# more info : http://www.sthda.com/english/wiki/wiki.php?id_contents=7851
# Data treatment
pca <- prcomp(circ_counts3, center = TRUE, scale = TRUE)

# General plots
scree_plot         <- fviz_screeplot(pca, ncp = 10)
graph_of_variables <- fviz_pca_var(pca, col.var = "contrib") +
  theme_minimal()
# Individual plots
dim1_dim2 <- fviz_pca_ind(pca, col.ind="cos2") +
  scale_color_gradient2(low="blue", mid="white",
                    high="red", midpoint=0.50)
dim1_dim3 <- fviz_pca_ind(pca, axes = c(1,3), col.ind="cos2") +
  scale_color_gradient2(low="blue", mid="white",
                    high="red", midpoint=0.50)
dim2_dim3 <- fviz_pca_ind(pca, axes = c(2,3), col.ind="cos2") +
  scale_color_gradient2(low="blue", mid="white",
                    high="red", midpoint=0.50)

################################################################################
# 3. Save ggplots
################################################################################

plots <- list(heatmap = heatmap, boxplot = boxplot, violinplot = violinplot,
              histogram = histogram, scree_plot = scree_plot,
              graph_of_variables = graph_of_variables, dim1_dim2 = dim1_dim2,
              dim1_dim3 = dim1_dim3, dim2_dim3 = dim2_dim3)

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
