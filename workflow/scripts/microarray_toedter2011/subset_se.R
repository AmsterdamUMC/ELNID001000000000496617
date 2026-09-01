#!/usr/bin/env Rscript
# Subset the general SummarizedExperiment to the samples needed for the comparison

# Paths
se_rds <- snakemake@input[["se_rds"]]
filter_params <- snakemake@params[["filter"]]
se_subset_rds <- snakemake@output[["se_subset_rds"]]
log_file <- snakemake@log[[1]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Libraries
library(tidyverse)
library(SummarizedExperiment)
library(limma)
library(stringr)

sessionInfo()

# Import full dds object
toedter_se <- readRDS(se_rds)

selected_samples <- colData(toedter_se) %>%
  data.frame() %>%
  dplyr::filter(Tissue %in% filter_params$Tissue,
                Response %in% filter_params$Response,
                Week %in% filter_params$Week,
                Treatment %in% filter_params$Treatment)

se_subset <- toedter_se[,colnames(toedter_se) %in% rownames(selected_samples)]

assay(se_subset, "logqnorm") <- normalizeBetweenArrays(assay(se_subset, 
                                                             i = "log"),
                                                       method = "quantile")

# Save
saveRDS(object = se_subset, file = se_subset_rds, compress = "gzip")

sink(type = "message")
sink(type = "output")
close(log_file)