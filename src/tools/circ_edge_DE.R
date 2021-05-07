#!bin/bash/R

################################################################################
# Differential expression analysis of circular RNAs. Script prepared for work with
# CIRIquant output data.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 25-03-2021
# Last modification : 06-05-2021
################################################################################

suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("statmod"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))

# !! Hacer una nota al usuario de poner el nombre en de las columnas con la primera
# en mayúscula, siempre y cuando ponga él sus propios archivos
check_lib <- function(lib){
  if(ncol(lib) == 5 ){
    lib      <- lib %>% rename(Group = group)
    metadata <- as.factor(lib['group'])
    return(metadata)
  }elif(ncol(lib) == 6){
    lib            <- lib %>% rename(Group = group)
    lib            <- lib %>% rename(Sex = sex)
    metadata       <- cbind(lib['group'], lib['sex'])
    metadata$group <- as.factor(metadata$group)
    metadata$sex   <- as.factor(metadata$sex)
    return(metadata)
  }else{
    stop(paste0("ERROR: Invalid number of columns in library information file."))
  }
}


option_list <- list(
  make_option(c("-d", "--design"), action = "store", type = "character", default = NULL,
                help = "experimental design", metavar = "design"),
  make_option(c("--lib"), type = "character", default = NULL,
              help = "File containing library information about samples:
              sample names, total reads, mapped reads, circular reads, Group and Sex", metavar = "FILE1"),
  make_option(c("-c", "--circ_counts"), default = NULL, help = "Numeric matrix of circular counts",
              metavar="FILE2" ),
  make_option(c("-l", "--linear_counts"), default = NULL, help = "Numeric matrix of linear counts",
              metavar = "FILE3" ),
  #make_option(c("-g", "--gene"), default=NULL, help="gene information file", metavar="FILE4"),
   # mirar si aquesta opció es possible incorporarla en funció de la info que hagi en gene_matrix_counts.csv
  make_option(c("-C", "--circ_info"), default = NULL, help = "Circular information file,
              containning circRNA annotation information.", metavar = "FILE5"),
  make_option(c("-p", "--pvalue"), type = "numeric", default = 0.05,
              help = "pvalue", metavar = "pval" ),
  make_option(c("-f", "--foldchange"),type = "numeric", default = 1.5, help = "foldchange",
              metavar = "FC" ),
  make_option("--outdir"), default = NULL, help = "output file")
  )

parser <- OptionParser(option_list=option_list)
opt    <- parse_args(parser)

if(is.null(opt$design)|| is.null(opt$lib) || is.null(opt$circ_counts) ||
  is.null(opt$linear_counts) || is.null(opt$circ_info)){
  print_help(parser)
  stop("Options --design/--lib/--circ_counts/--linear_counts/--circ_info
   must be supplied\n", call.=FALSE)

}else{
  cat("The supplied arguments are the followings:\n")
  cat(paste(" Experimental design      = ", opt$design, "\n",
            "Library file              = ", opt$lib, "\n",
            "Circular count file       = ", opt$circ_counts, "\n",
            "Linear count file         = ", opt$linear_counts, "\n",
            #"Gene information file     = ", opt$gene, "\n",
            "Circular information file = ", opt$circ_info, "\n"
          )
        )
}

# 1. Load data
print("Loading data...")

lib_mtx        <- read.csv(opt$lib, row.names = 1)
metadata       <- check_lib(lib_mtx)
circ_info      <- read.csv(opt$circ_info)
circrna_counts <-as.matrix(read.csv(opt$circ_counts))
linear_counts  <- as.matrix(read.csv(opt$linear_counts))
linear_counts  <- linear_counts[ , rownames(lib_mtx)]
#gene_info <- read.csv(opt$gene) # Annotation information for each gene. Segurament la obtindre de la matriu de linear counts
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
circrna_df     <- cbind(circ_order, circrna_df)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~VOLCANO PLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Data treatment: creating groups
pval    <- -log2(opt$pvalue)
FC      <- log2(opt$foldchange)

DE_data["minusLog2Pvalue"] <- -log2(DE_data$PValue)
DE_data["group"] <- "Not_Significant"
DE_data[which(DE_data['minusLog2Pvalue'] > pval & DE_data["logFC"] > FC ),"group"] <- "Up"
DE_data[which(DE_data['minusLog2Pvalue'] > pval & DE_data["logFC"] < -FC ),"group"] <- "Down"
DE_data$group <- as.factor(DE_data$group)

#Top20 circRNAs (by pvalue), which names will be ploted
top20   <- DE_data$id[1:20]
DE_data <- DE_data %>% mutate(plotname = ifelse(id %in% top20, name, "" ))

# Plot
volcano_plot <-
 ggplot(DE_data, aes(x = logFC, y = minusLog2Pvalue, color = group)) +
    geom_point(size = 0.5)+
    scale_colour_manual(values = c("red", "grey", "blue"))+
    geom_text_repel(aes(label = plotname), size = 3)+
    geom_hline(yintercept = pval, linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = c(FC, -FC), linetype = "dashed", size = 0.5 ) +
    xlab("log2FC") +
    ylab("-log10(p-value)")

# Download data
DE_matrix_path <- paste0(opt$outdir,"/circrna_DE.csv")
volcano_path   <- paste0(opt$outdir,"/volcano_plot.svg")

write.csv(circrna_df, file = DE_matrix_path, quote = FALSE)
ggsave(filename = volcano_path , plot = volcano_plot, device = "svg")

print(paste("Differential Expression analysis done. Output files:\n",
            DE_matrix_path, "\n", volcano_path, "\n"))
