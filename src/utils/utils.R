#!/bin/R

###############################################################################
# R script containning useful functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 10-05-2021
# Last modification : 11-05-2021
###############################################################################

# Dependencies
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("dplyr"))

# Functions
# !! Hacer una nota al usuario de poner el nombre en de las columnas con la primera
# en mayúscula, siempre y cuando ponga él sus propios archivos

check_metadata <- function(opt_metadata, sep){
  metadata    <- read.csv(opt_metadata, sep = sep)
  if('Sample' %in% colnames(metadata)){
    metadata          <- metadata %>% dplyr::rename(sample = Sample)
    rownames(metadata)<- metadata$sample
  }else{
    stop(print("ERROR: 'Sample' column must be provided in metadata."))
  }
  if('Group' %in% colnames(metadata) & 'Sex' %in% colnames(metadata) ){
    metadata       <- metadata %>% dplyr::rename(group = Group)
    metadata       <- metadata %>% dplyr::rename(sex = Sex)

    metadata$group <- as.factor(metadata$group)
    metadata$sex   <- as.factor(metadata$sex)
    return(metadata)
  }else if('Group' %in% colnames(metadata)){
    metadata       <- metadata %>% dplyr::rename(group = Group)
    metadata$group <- as.factor(metadata$group)
    return(metadata)
  }else{
    stop(paste0("ERROR: 'Group' column must be supplied in metadata."))
  }
}


check_norm <- function(norm, metadata, design){
  design <- as.formula(design)
  if(norm == "False"){
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

check_data_co <- function(circ_counts, circ_info){
  if('id' %in% colnames(circ_info) & 'start' %in% colnames(circ_info) & 'end' %in% colnames(circ_info) ){
    circ_counts           <- as.matrix(read.csv(circ_counts))
    rownames(circ_counts) <- circ_info$id
    return(circ_counts)
  }else{
    stop(print("ERROR: 'id'/'start'/'end' columns must be present in gene_length file."))
  }
}

check_sep <- function(separator){
  switch(separator,
    comma     = separator <- ',' ,
    semicolon = separator <- ';' ,
    tab       = separator <- '\t',
    whitespace= separator <- ' ' )
  return(separator)
}
