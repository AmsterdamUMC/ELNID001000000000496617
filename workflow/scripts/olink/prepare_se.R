#!/usr/bin/env Rscript
# Build the general SummarizedExperiments from the NPX values and sample metadata

# Paths
olink_files_xlsx <- snakemake@input[["olink_files_xlsx"]]
sample_metadata_xlsx <- snakemake@input[["sample_metadata_xlsx"]]
se_rds <- snakemake@output[["se_rds"]]
log_file <- snakemake@log[[1]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Libraries
library(tidyverse)
library(SummarizedExperiment)
library(readxl)

sessionInfo()

# Preparing npx
olink_files <- readxl::read_excel(olink_files_xlsx)
npx_files <- unique(olink_files$Basename)

npx_full <- do.call(rbind, lapply(npx_files, function(npx_file){
  npx <- read.csv(npx_file, sep = ",")
  return(npx)
})) %>%
  dplyr::rename(ClientAccessionID = SampleID) %>%
  dplyr::inner_join(olink_files, by = "ClientAccessionID")

# Prepare expression
message("--- Preparing expression ---")
npx <- npx_full %>%
  dplyr::select(OlinkID, SampleID, NPX) %>%
  tidyr::pivot_wider(id_cols = OlinkID, names_from = SampleID, values_from = NPX) %>%
  tibble::column_to_rownames(var = "OlinkID")

# Prepare feature metadata
message("--- Preparing feature metadata ---")
feature_metadata <- npx_full %>%
  dplyr::select(OlinkID, UniProt, Assay, Panel, Panel_Lot_Nr) %>%
  unique() %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(var = "OlinkID")

# Prepare sample metadata
message("--- Preparing sample metadata ---")
sample_metadata <- readxl::read_excel(sample_metadata_xlsx) %>%
  tibble::column_to_rownames(var = "SampleID") %>%
  dplyr::mutate(Week = factor(Week, levels = c("W0", "W4", "W8", "EW")),
                Visit = factor(Visit, levels = c("Baseline", "Week 4", "Week 8", "Early Withdrawal")),
                Response = factor(Response, levels = c("NR", "LR", "R")),
                Age = as.numeric(Age),
                wk0_eMayo = as.numeric(wk0_eMayo),
                wk8_eMayo = as.numeric(wk8_eMayo),
                wk0_totalMayo = as.numeric(wk0_totalMayo),
                wk8_totalMayo = as.numeric(wk8_totalMayo),
                wk0_RHI = as.numeric(wk0_RHI),
                wk8_RHI = as.numeric(wk8_RHI),
                wk0_crp = as.numeric(wk0_crp))

# Prepare SummarizedExperiment
message("--- Creating SummarizedExperiment ---")
olink_se <- SummarizedExperiment(assays = list(npx = npx),
                                 rowData = feature_metadata[rownames(npx),],
                                 colData = sample_metadata[colnames(npx),])

# Save
saveRDS(olink_se, se_rds, compress = "gzip")
message("Saved: ", se_rds)

sink(type = "output")
sink(type = "message")
close(log_file)