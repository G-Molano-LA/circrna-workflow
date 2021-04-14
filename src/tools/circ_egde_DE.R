#!bin/bash/R

###############################################################################
## Differential expression analysis
###############################################################################
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("statmod"))
suppressPackageStartupMessages(library("optparse"))

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
#               help="Show this help message and exit")

option_list <- list(
  make_option(c("-d", "--design"), action="store", type="character", default=NULL,
                help="experimental design", metavar="design"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="metadata file", metavar="FILE1" ),
  make_option(c("-c", "--circ_counts"), default=NULL, help="circular counts file",
              metavar="FILE2" ),
  make_option(c("-l", "--linear_counts"), default=NULL, help="linear counts file",
              metavar="FILE3" ),
  make_option(c("-g", "--gene"), default=NULL, help="gene information file",
              metavar="FILE4"),
  make_option(c("-C", "--circ_info"), default=NULL, help="circular information file",
              metavar="FILE5"),
  make_option(c("-o", "--out"), default="circrna_edge_DE.csv", help="output file",
              metavar="output")
  )

parser=OptionParser(option_list=option_list)
opt=parse_args(parser)

if(is.null(opt$design)|| is.null(opt$metadata) || is.null(opt$circ_counts) ||
  is.null(opt$linear_counts) || is.null(opt$gene) || is.null(opt$circ_info)){
  print_help(parser)
  stop("Options --design/--metadata/--circ_counts/--linear_counts/--gene/--circ_info
   must be supplied\n", call.=FALSE)

}else{
  cat("The supplied arguments are the followings:\n")
  cat(paste(" Experimental design      = ", opt$design, "\n",
            "Metadata file             = ", opt$metadata, "\n",
            "Circular count file       = ", opt$circ_counts, "\n",
            "Linear count file         = ", opt$linear_counts, "\n",
            "Gene information file     = ", opt$gene, "\n",
            "Circular information file = ", opt$circ_info, "\n",
            "Output file               = ", opt$out, "\n"
          )
        )
}

# 1. Load data
print("Loading data...")
metadata<-read.csv(opt$metadata,stringsAsFactors=TRUE) # data frame containing information
                                                      #+ for each sample
# Make sure that in your metadata information there is a column names' condition'

gene_info <- read.csv(opt$gene) # data frame containing annotation
                                            #+ information for each gene
circ_info <- read.csv(opt$circ_info) # data frame containing annotation
                                            #+ information for each circular RNA
circrna_counts <-as.matrix(read.csv(opt$circ_counts))  # numeric matrix of circular counts
linear_counts <- as.matrix(read.csv(opt$linear_counts)) # numeric matrix of linear counts
# shell: wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159225/suppl/GSE159225%5FcircRNA%5FReadcount%5FAllsamples%2Ecsv%2Egz
#         -O docs/circrna_counts.csv
#        wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159225/suppl/GSE159225%5FlinearRNA%5FReadcount%5FAllsamples%2Ecsv%2Egz
#         -O docs/linear_counts.csv

print("Done.")
# 2. The DGEList data class
## Linears
print("Normalizing circular counts by effective library sizes... The scale factor uses a trimmed mean of M-values (TMM) between each pair of samples.")

linear_DGE <- DGEList(counts   = linear_counts,
                      samples  = metadata,
                      genes    = gene_info)
linear_keep <- filterByExpr(linear_DGE)
linear_DGE <- linear_DGE[linear_keep, keep.lib.sizes=FALSE]
linear_DGE <- calcNormFactors(linear_DGE) # Normalizing by effective library
#+ sizes (gene expression level). The scale factor uses a trimmed mean of M-values
#+ (TMM) between each pair of samples --> Trimmed mean of log FC normalization

opt$design <- as.formula(opt$design)
design <- model.matrix (object=opt$design, data=linear_DGE$samples)
print("Done design")
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
print(paste("DE analysis done. Output file created: ", opt$out))
