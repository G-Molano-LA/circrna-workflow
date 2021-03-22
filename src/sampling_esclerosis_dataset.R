#!usr/bin/R

################################################################################
## R script for sampleling
## Author: G. Molano, LA
################################################################################

# 1. Metadata
metadata<-read.csv("~/Desktop/TFG/data/esclerosis/processed_data/metadata_SRA.txt",
                   stringsAsFactors=TRUE) # Download from SRA repository
metadata_filtered <- subset(metadata,  select=c(Run, disease_state, ms_type, sex))
head(metadata_filtered)

# 2. Sampling metadata
set.seed(108)
metadata_sample <- metadata_filtered %>% group_by(disease_state, sex) %>% slice_sample(n=2)
metadata_sample

# 3. Download metadata_sample
write.csv(metadata_sample, file="metadata_sampled.csv", row.names=F)
