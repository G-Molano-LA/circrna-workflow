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

# 1. Input
option_list <- list(
  make_option(c("-p", "--pvalue"), type="numeric", default=0.05,
              help="pvalue", metavar="pval" ),
  make_option(c("-f", "--foldchange"),type="numeric", default=1.5, help="foldchange",
              metavar="FC" ),
  make_option(c("-o", "--output"),type="character", default="pdf")
  )
parser=OptionParser(option_list=option_list)
opt=parse_args(parser)


# 2. Data treatment
DE_data <- read.csv("libs/DE_analysis/circrna_DE.csv", row.names=1)
pval = -log2(opt$pvalue)
FC = log2(opt$foldchange)

DE_data["minusLog2Pvalue"] <- -log2(DE_data$PValue)
DE_data["group"] <- "Not_Significant"
DE_data[which(DE_data['minusLog2Pvalue'] > pval & DE_data["logFC"] > FC ),"group"] <- "Up"
DE_data[which(DE_data['minusLog2Pvalue'] > pval & DE_data["logFC"] < -FC ),"group"] <- "Down"
DE_data$group <- as.factor(DE_data$group)

# 2.1. Top10 circRNAs (by pvalue), which names will be ploted
top20<- DE_data$id[1:20]
DE_data <- DE_data %>%
  mutate(plotname=ifelse(id %in% top20, name, "" ))

# 3. VOLCANO PLOT
volcano_plot = function(data){
  df <- DE_data
  plot <- ggplot(df, aes(x = logFC, y = minusLog2Pvalue, color = group)) +
    geom_point(size=0.5)+
    scale_colour_manual(values=c("red", "grey", "blue"))+
    geom_text_repel(aes(label= plotname), size=3)+
    geom_hline(yintercept = pval, linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = c(FC, -FC), linetype = "dashed", size = 0.5 ) +
    xlab("log2FC") +
    ylab("-log10(p-value)")

    switch(opt$output,
      pdf = ggsave("libs/plots/volcano_plot.pdf",plot),
      svg = ggsave("libs/plots/volcano_plot.svg",plot)
    )

}

volcano_plot(data=DE_data)

# 4. HEATMAP
# 4.1. Load NORMALIZED circular count matrix
norm_counts <- read.csv("libs/network/TMM_counts.txt", sep = "", row.names=1)

# 4.2. Plot top10 DE (by pvalue) circular RNAs
norm_counts1 <- filter(norm_counts, rownames(norm_counts) %in% top20)
norm_counts1 <- as.matrix(norm_counts1)
norm_counts1 <- melt(norm_counts1)
colnames(norm_counts1) <- c("Circular_RNAs", "Samples", "value")

ggplot(norm_counts1, aes(x=Samples, y=Circular_RNAs, fill=value)) +
  geom_tile(colour = "white") +
  labs(fill = "TMM counts") +
  scale_fill_gradient(low="white", high = "steelblue") +
  geom_text(aes(label = round(value, 3)))

# 5. BOXPLOT
# 5.1
ggplot(norm_counts1, aes(x = factor(Samples), y = value, fill = factor(Samples))) +
  geom_boxplot() +
  geom_jitter(color="grey", size=0.4, alpha=0.9) +
  labs(fill = "Samples") +
  xlab("Samples") +
  ylab("Normalized TMM counts")

# 6. VIOLIN PLOT
ggplot(norm_counts1, aes(x=factor(Samples), y = value, fill = factor(Samples))) +
  geom_violin() +
  labs(fill = "Samples") +
  xlab("Samples") +
  ylab("Normalized TMM counts")


ggplot(norm_counts1, aes(x=value)) +
  geom_histogram() +
  xlab("Number of reads")
