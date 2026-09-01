#!/usr/bin/env Rscript
# Import and create a pseudobulk counts. No further filtering is performed to mimick a bulk RNA-sequencing experiment.

# Paths
mtx_dir <- snakemake@input[["melonardanaz2025_dir"]]
pbcounts_csv <- snakemake@output[["pbcounts_csv"]]
sampleid <- snakemake@params[["sampleid"]]
log_file <- snakemake@log[[1]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Libraries
library(Seurat)

sessionInfo()

# Preparing counts
message("--- Preparing counts ---")
counts10x <- Read10X(mtx_dir)
pbcounts <- rowSums(counts10x)
pbcounts <- data.frame(pbcounts)
colnames(pbcounts) <- sampleid

# Save
message("--- Saving pseudobulk counts ---")
write.csv(pbcounts, pbcounts_csv, row.names = T)

sink(type = "message")
sink(type = "output")
close(log_file)