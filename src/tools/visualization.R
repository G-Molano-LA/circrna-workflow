#!/bin/R

# Dependencies
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("VennDiagram"))

# 1. Input
option_list <- list(
  make_option(c("-d", "--data"), default=NULL, help="data file", metavar="DATA"),
  make_option(c("-p", "--pvalue"), type="numeric", default=0.05,
              help="pvalue", metavar="pval" ),
  make_option(c("-f", "--foldchange"),type="numeric", default=1.5, help="foldchange",
              metavar="FC" ),
  make_option(c("-o", "--output"),type="character", default="pdf")
  )
parser=OptionParser(option_list=option_list)
opt=parse_args(parser)

if(is.null(opt$data)){
  print_help(parser)
  stop("Options --data must be supplied.\n", call.=FALSE)
}

# 2. Data treatment
data <- read.csv(opt$data, row.names=1)
pval = -log2(opt$pvalue)
FC = log2(opt$foldchange)

data["minusLog2Pvalue"] <- -log2(data$PValue)
data["group"] <- "Not_Significant"
data[which(data['minusLog2Pvalue'] > pval & data["logFC"] > FC ),"group"] <- "Up"
data[which(data['minusLog2Pvalue'] > pval & data["logFC"] < -FC ),"group"] <- "Down"
data$group <- as.factor(data$group)

top10<- data$logFC[1:10]
data <- data %>%
  mutate(plotname=ifelse(logFC %in% top10, name, "" ))


# 3. Volcano Plot

volcano_plot = function(data){
  df <- data
  plot <- ggplot(df, aes(x = logFC, y = minusLog2Pvalue, color = group)) +
    geom_point(size=0.5)+
    scale_colour_manual(values=c("red", "grey", "blue"))+
    geom_text_repel(aes(label= plotname), size=3)+
    geom_hline(yintercept = pval, linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = c(FC, -FC), linetype = "dashed", size = 0.5 ) +
    xlab("log2FC") +
    ylab("-log10(p-value)")

    switch(opt$output,
      pdf = ggsave("volcano_plot.pdf",plot),
      svg = ggsave("volcano_plot.svg",plot)
    )

}

volcano_plot(data=data)

# 4. Venn Diagram
# venn.diagram(
#   x=list(ciri_list, circexplorer-list), # data from ciri and circexplorer merged results
#   category.names=c("CIRI2", "CircExplorer2"),
#   filename="venn_diagram.png",
#   output=TRUE,
#   fill=c("lightblue", "palegreen")
# )

# 5. HeatMap
# Load data
# data <- read.csv("~/circrna-workflow/docs/circrna_counts.csv", sep=";")
#     # all counts: not available yet --> output from merge ciri and cirexplorer results
#     #           or from merged results from ciriquant
#
# top10 <- data$id[1:10] # data from DE analysis
# data2 <- filter(data, id %in% top10)
# colnames <- data2$id
# data2 <- select(data2, !c(id, X.chrom, start, end, score, strand, name))
#
# data2 <- as.matrix(data2)
# data2 <- t(data2)
# colnames(data2) <- colnames
#
# heatmap(data2, Colv = NA, Rowv = NA, scale = "column") # scaling by genes
