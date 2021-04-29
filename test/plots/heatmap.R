
library("dplyr")
library("ggplot2")
library("reshape2")

# R BASIC
# Load data : REVISAR NORMALIZACIÓN DE DATOS
circ_counts <- read.csv("libs/network/TMM_counts.txt", sep = "", row.names=1)
de_data <- read.csv("~/circrna-workflow/libs/DE_analysis/circrna_edge_DE.csv")
# Plotting the top 10 DE (by pvalue) circRNAs
top20 <- de_data$id[1:20]
circ_counts <- filter(circ_counts, rownames(circ_counts) %in% top20)
circ_counts1 <- as.matrix(circ_counts)


# Normalizing counts in each sample
#heatmap(circ_counts1, scale = "column")

# GGPLOT2
circ_counts2 <- melt(circ_counts1)
colnames(circ_counts2) <- c("Circular_RNAs", "Samples", "value")
ggplot(circ_counts2, aes(x=Samples, y=Circular_RNAs, fill=value)) +
  geom_tile(colour = "white") +
  labs(fill = "TMM counts") +
  scale_fill_gradient(low="white", high = "steelblue") +
  geom_text(aes(label = round(value, 3)))
