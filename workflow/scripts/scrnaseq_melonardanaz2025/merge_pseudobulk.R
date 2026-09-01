#!/usr/bin/env Rscript
# Merge the SeuratObjects into a single large SeuratObject

# Paths
pbcounts_csvs <- snakemake@input[["pbcounts_csvs"]]
pbcounts_merged_csv <- snakemake@output[["pbcounts_merged_csv"]]
log_file <- snakemake@log[[1]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Libraries
library(dplyr)

sessionInfo()

# print(pbcounts_csvs)

pbcounts <- lapply(pbcounts_csvs, function(pbcounts_csv){
  pbcount <- read.csv(pbcounts_csv) %>%
    tibble::column_to_rownames(var = "X")
  return(pbcount)
})
genes <- unique(unlist(lapply(pbcounts, rownames)))
print(genes)
sampleids <- unlist(lapply(pbcounts, colnames))
pbcounts_merged <- do.call(cbind, lapply(pbcounts, function(sampleid){
  sampleid[genes,]
}))
rownames(pbcounts_merged) <- genes
colnames(pbcounts_merged) <- sampleids

# Save
message("--- Saving merged pseudobulk counts ---")
write.csv(pbcounts_merged, pbcounts_merged_csv)

sink(type = "message")
sink(type = "output")
close(log_file)