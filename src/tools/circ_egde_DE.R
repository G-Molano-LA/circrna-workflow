#!bin/bash/R

###############################################################################
## Differential expression analysis
###############################################################################
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("statmod"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
#               help="Show this help message and exit")

option_list <- list(
  make_option(c("-d", "--design"), action="store", type="character", default=NULL,
                help="experimental design", metavar="design"),
  make_option(c("--lib"), type="character", default=NULL,
              help="File containing library information about samples:
              sample names, total reads, mapped reads, circular reads, Group and Sex", metavar="FILE1" ),
  make_option(c("-c", "--circ_counts"), default=NULL, help="Numeric matrix of circular counts",
              metavar="FILE2" ),
  make_option(c("-l", "--linear_counts"), default=NULL, help="Numeric matrix of linear counts",
              metavar="FILE3" ),
  #make_option(c("-g", "--gene"), default=NULL, help="gene information file", metavar="FILE4"),
   # mirar si aquesta opció es possible incorporarla en funció de la info que hagi en gene_matrix_counts.csv
  make_option(c("-C", "--circ_info"), default=NULL, help="Circular information file,
    containning circRNA annotation information.",
              metavar="FILE5"),
  make_option(c("-o", "--out"), default="circrna_edge_DE.csv", help="output file",
              metavar="output")
  )

parser=OptionParser(option_list=option_list)
opt=parse_args(parser)

if(is.null(opt$design)|| is.null(opt$lib) || is.null(opt$circ_counts) ||
  is.null(opt$linear_counts) || is.null(opt$circ_info)){
  print_help(parser)
  stop("Options --design/--lib/--circ_counts/--linear_counts/--circ_info
   must be supplied\n", call.=FALSE)

}else{
  cat("The supplied arguments are the followings:\n")
  cat(paste(" Experimental design      = ", opt$design, "\n",
            "Library file             = ", opt$lib, "\n",
            "Circular count file       = ", opt$circ_counts, "\n",
            "Linear count file         = ", opt$linear_counts, "\n",
            #"Gene information file     = ", opt$gene, "\n",
            "Circular information file = ", opt$circ_info, "\n",
            "Output file               = ", opt$out, "\n"
          )
        )
}

# 1. Load data
print("Loading data...")

lib_mtx <-read.csv(opt$lib, row.names = 1)
#gene_info <- read.csv(opt$gene) # Annotation information for each gene. Segurament la obtindre de la matriu de linear counts
circ_info <- read.csv(opt$circ_info)
circrna_counts <-as.matrix(read.csv(opt$circ_counts))
linear_counts <- as.matrix(read.csv(opt$linear_counts))
linear_counts <- linear_counts[ , rownames(lib_mtx)]

lib_mtx<- lib_mtx %>% rename(Group = group)
lib_mtx<- lib_mtx %>% rename(Sex = sex)
metadata <- cbind(lib_mtx['group'], lib_mtx['sex'])
metadata$group <- as.factor(metadata$group)
metadata$sex <- as.factor(metadata$sex)

print("Done.")

# 2. The DGEList data class
# 2.1. Linears
print("Normalizing circular counts by effective library sizes... The scale factor uses a trimmed mean of M-values (TMM) between each pair of samples.")

linear_DGE <- DGEList(counts   = linear_counts,
                      samples = metadata)
                      #genes    = gene_info)
linear_keep <- filterByExpr(linear_DGE)
linear_DGE <- linear_DGE[linear_keep, keep.lib.sizes=FALSE]
linear_DGE <- calcNormFactors(linear_DGE) # Normalizing by effective library
#+ sizes (gene expression level). The scale factor uses a trimmed mean of M-values
#+ (TMM) between each pair of samples --> Trimmed mean of log FC normalization

opt$design <- as.formula(opt$design)
design <- model.matrix (object=opt$design, data=linear_DGE$samples)

## Circulars
circrna_DGE <- DGEList(counts = circrna_counts,
                       samples = metadata,
                       genes= circ_info,
                       lib.size = linear_DGE$samples[, "lib.size"],
                       norm.factors = linear_DGE$samples[, "norm.factors"])

print("Done.")
print("Starting differential expression analysis of circular RNAs...")

circrna_DGE <- estimateDisp(circrna_DGE, design, robust = TRUE)
circrna_fit <- glmFit(circrna_DGE, design)
circrna_lrt <- glmLRT(circrna_fit)

print("Summary:")
summary(decideTests(circrna_lrt))


circrna_df <- circrna_lrt$table
pval_order <- order(circrna_lrt$table$PValue)
circrna_df$DE <- as.vector(decideTestsDGE(circrna_lrt))
circrna_df <- circrna_df[pval_order, ]
circrna_df$FDR <- p.adjust(circrna_df$PValue, method="fdr")
circ_order <- circrna_lrt$genes[pval_order, ]
circrna_df <-cbind(circ_order, circrna_df)
write.csv(circrna_df, file =opt$out, quote=FALSE)
print(paste("Differential Expression analysis done. Output file created: ", opt$out))
