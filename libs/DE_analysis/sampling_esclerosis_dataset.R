#!usr/bin/R

################################################################################
## R script for sampleling
## Author: G. Molano, LA
################################################################################
library(dplyr)

# 1. Metadata
metadata<-read.csv("~/alejandra/circrna-workflow/docs/metadata_all.txt",
                   stringsAsFactors=TRUE) # Download from SRA repository
metadata_filtered <- subset(metadata,  select=c(Run, disease_state, ms_type, sex))
metadata_filtered <- metadata_filtered %>% rename(group = disease_state) # must be a colunm name group for edgeR analysis
metadata_filtered <- metadata_filtered[order(metadata_filtered$Run),]

## 2. Adding patient id
rrms <-character()
for(i in 1:20){
  if(i<10){
    rrms[i]<-paste0("RRMS0",i)
  }else{
    rrms[i]<-paste0("RRMS",i)
  }
}

spms <-character()
for(i in 1:10){
  if(i==10){
    spms[i]<-paste0("SPMS",i)
  }else{
    spms[i]<-paste0("SPMS0",i)
  }
}

hc <-character()
for(i in 1:20){
  if(i<10){
    hc[i]<-paste0("HC0",i)
  }else{
    hc[i]<-paste0("HC",i)
  }
}

patient_id <- c(rrms, spms, hc)
print(patient_id)
metadata_filtered$patient <- patient_id

# 2. Sampling metadata
set.seed(100)
metadata_sample <- metadata_filtered %>% group_by(group, sex) %>% slice_sample(n=2)
metadata_sample

# 3. Download metadata_sample
write.csv(metadata_sample, file="~/alejandra/circrna-workflow/docs/metadata_samples.csv", row.names=F)
