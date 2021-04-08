

source = "https://github.com/czllab/NetMiner/blob/master/manual.txt"

# Requirements
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")

bioc_pkg <- c("impute", "preprocessCore", "GO.db")
BiocManager::install(bioc_pkg)

pkg <- c("WGCNA", "corpcor", "bc3net", "igraph")
install.packages(pkg, dependencies = TRUE)
