#!/bin/R

###############################################################################
## R script for visualization data
## Author: G. Molano, LA
###############################################################################

# Dependencies
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("factoextra"))

################################################################################
# 1. Input
################################################################################
option_list <- list(
  make_option(c("-p", "--pvalue"), type="numeric", default=0.05,
              help="pvalue", metavar="pval" ),
  make_option(c("-f", "--foldchange"),type="numeric", default=1.5, help="foldchange",
              metavar="FC" ),
  make_option(c("-o", "--output"),type="character", default="pdf")
  )
parser=OptionParser(option_list=option_list)
opt=parse_args(parser)

################################################################################
# 2. Data: differential expression matrix + input values
################################################################################
DE_data <- read.csv("libs/DE_analysis/circrna_DE.csv", row.names=1)
pval = -log2(opt$pvalue)
FC = log2(opt$foldchange)

#~~~~~~~~~~~~~~~~~~~~~VOLCANO PLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Data treatment: creating groups
DE_data["minusLog2Pvalue"] <- -log2(DE_data$PValue)
DE_data["group"] <- "Not_Significant"
DE_data[which(DE_data['minusLog2Pvalue'] > pval & DE_data["logFC"] > FC ),"group"] <- "Up"
DE_data[which(DE_data['minusLog2Pvalue'] > pval & DE_data["logFC"] < -FC ),"group"] <- "Down"
DE_data$group <- as.factor(DE_data$group)
##Top10 circRNAs (by pvalue), which names will be ploted
top20<- DE_data$id[1:20]
DE_data <- DE_data %>%
  mutate(plotname=ifelse(id %in% top20, name, "" ))

# Plot
volcano_plot <-
 ggplot(DE_data, aes(x = logFC, y = minusLog2Pvalue, color = group)) +
    geom_point(size=0.5)+
    scale_colour_manual(values=c("red", "grey", "blue"))+
    geom_text_repel(aes(label= plotname), size=3)+
    geom_hline(yintercept = pval, linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = c(FC, -FC), linetype = "dashed", size = 0.5 ) +
    xlab("log2FC") +
    ylab("-log10(p-value)")


################################################################################
# 3. Data: normalized expression matrix
################################################################################
circ_counts <- read.csv("libs/network/TMM_counts.txt", sep = "", row.names=1)

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
dist1 <- dist(circ_counts3, method = "euclidean") # Euclidean distance
dendro <- hclust(dist1, method = "average") #  Hierarchical Clustering with hclust

# Save Plot
dendro <- as.dendrogram(dendro)
upper_ylim <- round(max(dist1), digits = -2) # round to hundredths

svg(file="/libs/plots/dendrograma.svg")
plot(dendro, ylab = "Euclidean Distance",xlab= "Samples", main = "Cluster Dendrogram",
  ylim = c(0,upper_ylim))
dev.off()

#~~~~~~~~~~~~~~~~~~~~~PRINCIPAL COMPONENTS ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# more info : http://www.sthda.com/english/wiki/wiki.php?id_contents=7851
# Data treatment
pca <- prcomp(circ_counts3, center = TRUE, scale = TRUE)

# General plots
scree_plot <- fviz_screeplot(pca, ncp=10)
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
# 4. Save ggplots
################################################################################

plots <- list(volcano_plot = volcano_plot,heatmap = heatmap, boxplot = boxplot,
              violinplot = violinplot, histogram = histogram, scree_plot = scree_plot,
              graph_of_variables = graph_of_variables, dim1_dim2 = dim1_dim2,
              dim1_dim3 = dim1_dim3, dim2_dim3 = dim2_dim3) # falta dendro

for (i in 1:length(plots)){
  switch(opt$output,
    svg = ggsave(filename = paste0("libs/plots/", names(plots[i]), ".svg"),
                 plot = plots[[i]], device = "svg"),
    pdf = ggsave(filename = paste0("libs/plots/", names(plots[i]), ".pdf"),
                 plot = plots[[i]], device = "pdf")
  )
}

switch(opt$output,
       svg = ggsave(filename = paste0("libs/plots/", names(plots[i]), ".svg"),
                    plot = plots[[i]], device = "svg"),
       pdf = ggsave(filename = paste0("libs/plots/", names(plots[i]), "pdf"),
                    plot = plots[[i]], device = "pdf")
)
