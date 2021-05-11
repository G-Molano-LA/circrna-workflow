#!/bin/R

###############################################################################
# R script containning useful functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: G. Molano, LA (gonmola@hotmail.es)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Date              : 10-05-2021
# Last modification : 10-05-2021
###############################################################################

# Dependencies
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("dplyr"))

# Functions
# !! Hacer una nota al usuario de poner el nombre en de las columnas con la primera
# en mayúscula, siempre y cuando ponga él sus propios archivos

check_metadata <- function(metadata){
  if('Sample' %in% colnames(metadata)){
    metadata        <- metadata %>% rename(sample = Sample)
  }else{
    stop(print("ERROR: 'Sample' column must be provided in metadata."))
  }
  if('Group' %in% colnames(metadata) & 'Sex' %in% colnames(metadata) ){
    metadata       <- metadata %>% rename(group = Group)
    metadata       <- metadata %>% rename(sex = Sex)

    metadata$group <- as.factor(metadata$group)
    metadata$sex   <- as.factor(metadata$sex)
    return(metadata)
  }else if('Group' %in% colnames(metadata)){
    metadata       <- metadata %>% rename(group = Group)
    metadata       <- as.factor(metadata['group'])
    return(metadata)
  }else{
    stop(paste0("ERROR: 'Group' column must be supplied in metadata."))
  }
}


check_norm <- function(norm, metadata, circ_info, design){
  design <- as.formula(design)

  if(norm == "False" && circ_info != 'None'){
    circ_info             <- read.csv(opt$circ_info, row.names = 1)

    rownames(circ_counts) <- rownames(circ_info)
    colnames(circ_counts) <- rownames(metadata)
    DESeq_count_data      <- DESeqDataSetFromMatrix(countData = circ_counts,
                                                    colData   = metadata,
                                                    design    = design)
    DESeq_count_data <- estimateSizeFactors(DESeq_count_data)
    circ_counts      <- as.matrix(counts(DESeq_count_data, normalized = TRUE))
    return(circ_counts)
  }else if(norm == "False" && circ_info == 'None'){
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

check_sep <- function(separator){
  switch(separator,
    comma     = separator <- ',' ,
    semicolon = separator <- ';' ,
    tab       = separator <- '\t',
    whitespace= separator <- ' ' )
  return(separator)
}
