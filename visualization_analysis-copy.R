# Dependencies
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("factoextra"))

source("~/circrna-workflow/workflow/utils/utils.R")


# Load data
sep         <- check_sep("comma")
metadata    <- check_metadata("~/circrna-workflow/libs/DE_analysis/metadata.csv", sep)
circ_counts <- check_norm("~/circrna-workflow/results/circ_counts.txt", "False", metadata, "~group+sex")
#sub_set
# all_counts  <- circ_counts
# circ_counts <- circ_counts[, c(2,8,14,20,26,32,38,44,50)]
################################################################################
# 2. Visualization
################################################################################
circ_counts2 <- melt(circ_counts)
colnames(circ_counts2) <- c("Circular_RNAs", "Samples", "value")
#~~~~~~~~~~~~~~~~~~~~~BOXPLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot
boxplot <-
  ggplot(circ_counts2, aes(x = factor(Samples), y = value+0.1, fill = factor(Samples))) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  labs(fill = "Samples") +
  xlab("Samples") +
  ylab("Normalized counts")

#~~~~~~~~~~~~~~~~~~~~~VIOLINPLOT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
violinplot <-
  ggplot(circ_counts2, aes(x=factor(Samples), y = value+0.1, fill = factor(Samples))) +
  geom_violin() +
  scale_y_continuous(trans='log10', breaks = c(0.1,1.1,10.1,100.1), labels = c(0,1,10,100)) +
  labs(fill = "Samples") +
  xlab("Samples") +
  ylab("Normalized counts")

#~~~~~~~~~~~~~~~~~~~~~HISTOGRAM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
tab       <- prop(table(circ_counts2$value)) 
# maxval    <- round(max(tab),digits = -4)
# maxval_2  <- round(max(tab[tab != max(tab)]),digits = -3) 


histogram <-
  ggplot(circ_counts2, aes(x = value)) +
  geom_histogram(bins = 100, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  scale_y_continuous(trans = "log10") +
                     #breaks = c(1, 10, 100, maxval_2, maxval,labels = c(1, 10, 100, maxval_2, maxval)), 
  ylab("Relative frequency of circRNAs") +
  xlab("Number of reads")
