#!/usr/bin/env Rscript
# Subset the general DESeqDataSet to the data needed for the comparison

# Paths
dds_rds <- snakemake@input[["dds_rds"]]
filter_params <- snakemake@params[["filter"]]
dds_subset_rds <- snakemake@output[["dds_subset_rds"]]
vst_subset_rds <- snakemake@output[["vst_subset_rds"]]
log_file <- snakemake@log[[1]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Libraries
library(tidyverse)
library(DESeq2)
library(stringr)

sessionInfo()

# Import full dds object
dds <- readRDS(dds_rds)

selected_samples <- colData(dds) %>%
  data.frame() %>%
  dplyr::filter(Tissue %in% filter_params$Tissue,
                Response %in% filter_params$Response,
                Week %in% filter_params$Week)

dds_subset <- dds[,colnames(dds) %in% rownames(selected_samples)]
vst_subset <- vst(dds_subset)

# Save
saveRDS(object = dds_subset, file = dds_subset_rds, compress = "gzip")
saveRDS(object = vst_subset, file = vst_subset_rds, compress = "gzip")

sink(type = "message")
sink(type = "output")
close(log_file)