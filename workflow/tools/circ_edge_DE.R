#!bin/bash/R

################################################################################
# Differential expression analysis of circular RNAs. Script prepared for work with
# CIRIquant output data.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Last modification : 30-06-2021
################################################################################

suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("statmod"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("reshape2"))

source("utils/utils.R")

parser <- ArgumentParser(description = 'Differential Expression Analysis')
parser$add_argument("--design", action = "store", type = "character", default = NULL,
              help = "experimental design", metavar = "design")
parser$add_argument("--metadata", type = "character", default = NULL,
            help = "File containing library information about samples:sample names, total reads, mapped reads, circular reads, Group and Sex",
            metavar = "FILE1")
parser$add_argument("--sep", default = ",", action = 'store', help = "Separator")
parser$add_argument("--circ_counts", default = NULL, help = "Numeric matrix of circular counts",
            metavar="FILE2")
parser$add_argument("--linear_counts", default = NULL, help = "Numeric matrix of linear counts",
            metavar = "FILE3")
parser$add_argument("--circ_info", default = NULL, help = "Circular information file, containning circRNA annotation information.",
            metavar = "FILE5")
parser$add_argument("--pval", default = 0.05, help = "pvalue", metavar = "pval" )
parser$add_argument("--fc", default = 1.5, help = "foldchange", metavar = "FC" )
parser$add_argument("--outdir", default = NULL, help = "output file")

opt    <- parser$parse_args()

if(is.null(opt$design)|| is.null(opt$metadata) || is.null(opt$circ_counts) ||
  is.null(opt$linear_counts) || is.null(opt$circ_info)){
  parser$print_help()
  stop("Options --design/--metadata/--circ_counts/--linear_counts/--circ_info
   must be supplied\n", call.=FALSE)
}else{
  cat("The supplied arguments are the followings:\n")
  cat(paste(" Experimental design      = ", opt$design, "\n",
            "Library file              = ", opt$metadata, "\n",
            "Circular count file       = ", opt$circ_counts, "\n",
            "Linear count file         = ", opt$linear_counts, "\n",
            #"Gene information file     = ", opt$gene, "\n",
            "Circular information file = ", opt$circ_info, "\n"
          )
        )
}

# 1. Load data
print("Loading data...")
sep            <- check_sep(opt$sep)
metadata       <- check_metadata(opt$metadata, sep)
circ_info      <- read.csv(opt$circ_info)
circrna_counts <- as.matrix(read.csv(opt$circ_counts, row.names = 1))
linear_counts  <- as.matrix(read.csv(opt$linear_counts, row.names = 1))
print("Done.")

#~~~~~~~~~~~~~~~~~~~~~~~~~DIFFERENTIAL EXPRESSION ANALYSIS~~~~~~~~~~~~~~~~~~~~~~
# 2. The DGEList data class
# 2.1. Linears
print("Normalizing circular counts by effective library sizes...
  The scale factor uses a trimmed mean of M-values (TMM) between each pair of samples.")

linear_DGE <- DGEList(counts   = linear_counts,
                      samples = metadata)
                      #genes    = gene_info)
linear_keep <- filterByExpr(linear_DGE)
linear_DGE  <- linear_DGE[linear_keep, keep.lib.sizes = FALSE]
linear_DGE  <- calcNormFactors(linear_DGE) # Normalizing by effective library
#+ sizes (gene expression level). The scale factor uses a trimmed mean of M-values
#+ (TMM) between each pair of samples --> Trimmed mean of log FC normalization

opt$design <- as.formula(opt$design)
design     <- model.matrix (object = opt$design, data = linear_DGE$samples)

## Circulars
circrna_DGE <- DGEList(counts       = circrna_counts,
                       samples      = metadata,
                       genes        = circ_info,
                       lib.size     = linear_DGE$samples[, "lib.size"],
                       norm.factors = linear_DGE$samples[, "norm.factors"])

print("Done.")
print("Starting differential expression analysis of circular RNAs...")

circrna_DGE <- estimateDisp(circrna_DGE, design, robust = TRUE)
circrna_fit <- glmFit(circrna_DGE, design)
circrna_lrt <- glmLRT(circrna_fit)

print("Summary:")
summary(decideTests(circrna_lrt))


circrna_df     <- circrna_lrt$table
pval_order     <- order(circrna_lrt$table$PValue)
circrna_df$DE  <- as.vector(decideTestsDGE(circrna_lrt))
circrna_df     <- circrna_df[pval_order, ]
circrna_df$FDR <- p.adjust(circrna_df$PValue, method="fdr")

circ_order     <- circrna_lrt$genes[pval_order, ]
DE_data        <- cbind(circ_order, circrna_df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~VOLCANO PLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Data treatment: creating groups
pval    <- -log2(as.numeric(opt$pval))
FC      <- log2(as.numeric(opt$fc))

DE_data["minusLog2Pvalue"] <- -log2(DE_data$PValue)
DE_data["group"] <- "Not_Significant"
DE_data[which(DE_data['minusLog2Pvalue'] > pval & DE_data["logFC"] > FC ),"group"] <- "Up"
DE_data[which(DE_data['minusLog2Pvalue'] > pval & DE_data["logFC"] < -FC ),"group"] <- "Down"
DE_data$group <- as.factor(DE_data$group)

#Top20 circRNAs (by pvalue), which names will be ploted
top20   <- DE_data$circ_id[1:20]
DE_data <- DE_data %>% dplyr::mutate(plotname = ifelse(circ_id %in% top20, gene_id, "" ))

# Plot
volcano_plot <-
 ggplot(DE_data, aes(x = logFC, y = minusLog2Pvalue, color = group)) +
    geom_point(size = 0.5)+
    scale_colour_manual(values = c("red", "black", "blue"))+
    geom_text_repel(aes(label = plotname), size = 3)+
    geom_hline(yintercept = pval, linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = c(FC, -FC), linetype = "dashed", size = 0.5 ) +
    xlab("log2FC") +
    ylab("-log10(p-value)")

#~~~~~~~~~~~~~~~~~~~~~HEATMAP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Data treatment: Plot top20 DE (by pvalue) circular RNAs
circrna_counts <- as.data.frame(circrna_counts)
circ_counts1   <- as.matrix(filter(circrna_counts, rownames(circrna_counts) %in% top20))

# Download data
DE_matrix_path <- paste0(opt$outdir,"/circrna_DE.csv")
volcano_path   <- paste0(opt$outdir,"/volcano_plot.svg")
heatmap_path   <- paste0(opt$outdir,"/heatmap.svg")

write.csv(DE_data, file = DE_matrix_path, quote = FALSE)

ggsave(filename = volcano_path , plot = volcano_plot, device = "svg")

svg(file = heatmap_path)
heatmap(circ_counts1, scale = "none")
dev.off()

print(paste("Differential Expression analysis done. Output files:\n",
            DE_matrix_path, "\n", volcano_path, "\n", heatmap_path, "\n"))
