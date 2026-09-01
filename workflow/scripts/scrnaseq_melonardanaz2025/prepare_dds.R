#!/usr/bin/env Rscript
# Build the general DESeqDataSet from counts and sample metadata

# Paths
counts_csv <- snakemake@input[["counts_csv"]]
sample_metadata_xlsx <- snakemake@input[["sample_metadata_xlsx"]]
dds_rds <- snakemake@output[["dds_rds"]]
log_file <- snakemake@log[[1]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Libraries
library(tidyverse)
library(DESeq2)
library(readxl)

sessionInfo()

# Preparing counts
message("--- Preparing counts ---")
raw_counts <- read.csv(counts_csv, header = TRUE)

processed_counts <- raw_counts %>% 
  tibble::column_to_rownames(var = "X")
processed_counts <- processed_counts[rowSums(processed_counts) != 0, ]

message("Genes with non-zero counts: ", nrow(processed_counts))

# Preparing sample metadata
message("--- Preparing metadata ---")
sample_metadata <- read_excel(sample_metadata_xlsx) %>%
  dplyr::mutate(Sample_ID = gsub("-", ".", Sample_ID)) %>%
  dplyr::filter(Sample_ID %in% colnames(processed_counts)) %>%
  tibble::column_to_rownames(var = "Sample_ID")
sample_metadata <- sample_metadata[colnames(processed_counts),]

# Prepare DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = processed_counts,
  colData = sample_metadata,
  design = ~1
)

# Save
saveRDS(dds, dds_rds, compress = "gzip")
message("Saved: ", dds_rds)

sink(type = "message")
sink(type = "output")
close(log_file)
