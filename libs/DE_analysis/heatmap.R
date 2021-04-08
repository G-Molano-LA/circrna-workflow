
library(dplyr)

# Load data
data <- read.csv("~/circrna-workflow/docs/circrna_counts.csv", sep=";")
de_data <- read.csv("~/circrna-workflow/libs/DE_analysis/circrna_DE.csv")
head(de_data)

top10 <- de_data$name[1:10]
data2 <- filter(data, name %in% top10)
colnames <- data2$name
data2 <- select(data2, !c(id, X.chrom, start, end, score, strand, name)) 

data3 <- as.matrix(data2)
data3 <- t(data3)
colnames(data3) <- colnames

head(data3)
heatmap(data3, Colv = NA, Rowv = NA, scale = "column")
