#!/bin/R

###############################################################################
# R script for visualization data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 09-04-2021
# Last modification : 11-05-2021
###############################################################################

# Dependencies
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("factoextra"))

source("utils/utils.R")

################################################################################
# 1. Input
################################################################################
parser <- ArgumentParser(description = 'Visualization')
parser$add_argument( "--data", default = NULL, help = "Normalized Count Matrix")
parser$add_argument("--norm", default = "no", help = "Boolean variable")
parser$add_argument("--design", default = "~group", help = "Experimental design")
parser$add_argument("--metadata", default = NULL, help =
  "File containing library information about samples: sample names, total reads, mapped reads, circular reads, Group and Sex")
parser$add_argument("--sep", default = ",", action = 'store', help = "Separator")
parser$add_argument( "--output",type = "character", default = "svg")
parser$add_argument( "--units",type = "character", default = "in")
parser$add_argument( "--width",type = "character", default = "7")
parser$add_argument( "--height",type = "character", default = "7")
parser$add_argument("--outdir", type = "character", default = NULL)

opt <- parser$parse_args()

# Load data
sep         <- check_sep(opt$sep)
metadata    <- check_metadata(opt$metadata, sep)

circ_counts <- check_norm(opt$data, opt$norm, metadata, opt$design)

################################################################################
# 2. Visualization
################################################################################
circ_counts2 <- melt(circ_counts)
colnames(circ_counts2) <- c("Circular_RNAs", "Samples", "value")
#~~~~~~~~~~~~~~~~~~~~~BOXPLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot
boxplot <-
  ggplot(circ_counts2, aes(x = factor(Samples), y = value+0.1, fill = factor(Samples))) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  labs(fill = "Samples") +
  xlab("Samples") +
  ylab("Normalized counts")

#~~~~~~~~~~~~~~~~~~~~~VIOLINPLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
violinplot <-
  ggplot(circ_counts2, aes(x=factor(Samples), y = value+0.1, fill = factor(Samples))) +
  geom_violin() +
  scale_y_continuous(trans='log10', breaks = c(0.1,1.1,10.1,100.1), labels = c(0,1,10,100)) +
  labs(fill = "Samples") +
  xlab("Samples") +
  ylab("Normalized counts")

#~~~~~~~~~~~~~~~~~~~~~HISTOGRAM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
histogram <-
  ggplot(circ_counts2, aes(x = value)) +
  geom_histogram(aes(y = stat(count) / sum(count)), breaks = seq(0,10,0.5), fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  scale_y_continuous(trans = "identity", breaks = seq(0,1,0.1)) +
  xlim(NA, 10) +
  ylab("Relative frequency of circRNAs") +
  xlab("Number of reads") + 
  ggtitle("Circular Counts")

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

plots <- list(boxplot = boxplot, violin_plot = violinplot,
              histogram = histogram, scree_plot = scree_plot,
              graph_of_variables = graph_of_variables, dim1_dim2 = dim1_dim2,
              dim1_dim3 = dim1_dim3, dim2_dim3 = dim2_dim3)

for (i in 1:length(plots)){
  switch(opt$output,
    svg = ggsave(filename = paste0(opt$outdir, "/", names(plots[i]), ".svg"),
                 plot = plots[[i]], device = "svg", width = as.integer(opt$width),
                height = as.integer(opt$height), units = opt$units),
    pdf = ggsave(filename = paste0(opt$outdir, "/", names(plots[i]), ".pdf"),
                 plot = plots[[i]], device = "pdf", width = as.integer(opt$width),
                height = as.integer(opt$height), units = opt$units),
    jpeg = ggsave(filename = paste0(opt$outdir, "/", names(plots[i]), ".jpeg"),
                plot = plots[[i]], device = "jpeg", width = as.integer(opt$width),
               height = as.integer(opt$height), units = opt$units),
    pgn = ggsave(filename = paste0(opt$outdir, "/", names(plots[i]), ".pgn"),
                plot = plots[[i]], device = "png", width = as.integer(opt$width),
               height = as.integer(opt$height), units = opt$units)
  )
}
