

raw_counts <- read.csv("test/data/normalized_counts/Raw_Count.txt", sep = "\t")
median_counts <- read.csv("test/data/normalized_counts/Median_Count.txt", sep = "\t")
fpkm_counts <- read.csv("test/data/normalized_counts/FPKM_Count.txt", sep = "\t")
vst_counts <- read.csv("test/data/normalized_counts/VST_Count.txt", sep = "\t")
tmm_counts <- read.csv("test/data/normalized_counts/TMM_Count.txt", sep = "\t")
uq_counts <- read.csv("test/data/normalized_counts/UQ_Count.txt", sep = "\t")

files <- list(Raw = raw_counts, FPKM = fpkm_counts, Median = median_counts,
              VST = vst_counts, TMM = tmm_counts, UQ = uq_counts)


for(i in 1:length(files)){
  files[[i]] <- files[[i]][c(1:5000),]
  write.table(files[[i]],
              file = paste0("test/data/normalized_counts/subset/", names(files[i]), "_Count.txt"),
              row.names = TRUE, quote = FALSE, sep = "\t")
}
