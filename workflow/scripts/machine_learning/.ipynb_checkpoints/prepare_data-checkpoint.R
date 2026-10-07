#!/usr/bin/env Rscript
# Export ML feature matrices + labels from the prepared RNA-seq (dds) and Olink (se) objects, for pretreatment (baseline) response prediction.

# Paths
dds_rds <- snakemake@input[["dds_rds"]]
se_rds <- snakemake@input[["se_rds"]]
mlp <- snakemake@params[["ml_params"]]
labels_csv <- snakemake@output[["labels_csv"]]
tx_binf_rvnr_w0_data_csv <- snakemake@output[["tx_binf_rvnr_w0_data_csv"]]
tx_pb_rvnr_w0_data_csv <- snakemake@output[["tx_pb_rvnr_w0_data_csv"]]
px_binf_rvnr_w0_data_csv <- snakemake@output[["px_binf_rvnr_w0_data_csv"]]
px_pb_rvnr_w0_data_csv <- snakemake@output[["px_pb_rvnr_w0_data_csv"]]

log_file  <- snakemake@log[[1]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Libraries
library(tidyverse)
library(DESeq2)
library(SummarizedExperiment)
library(data.table)

sessionInfo()

timepoint_col <- mlp$week_column
timepoint <- mlp$week
label_col <- mlp$label_column
classes <- mlp$classes
pos_class <- mlp$positive_class

# Import
dds <- readRDS(dds_rds)
se <- readRDS(se_rds)

# TX
sample_metadata <- data.frame(colData(dds))
vsd <- DESeq2::vst(dds)
gex <- assay(vsd)

## TX_BINF_RvNR_W0 
TX_BINF_RvNR_W0_sample_metadata <- sample_metadata %>%
  dplyr::filter(!!sym(label_col) %in% classes,
                !!sym(timepoint_col) %in% timepoint)

TX_BINF_RvNR_W0_gex <- t(gex[,rownames(TX_BINF_RvNR_W0_sample_metadata)])

## TX_PB_RvNR_W0
TX_PB_RvNR_W0_sample_metadata <- sample_metadata %>%
  dplyr::filter(!!sym(label_col) %in% classes,
                !!sym(timepoint_col) %in% timepoint)

TX_PB_RvNR_W0_gex <- t(gex[,rownames(TX_PB_RvNR_W0_sample_metadata)])

# Olink
sample_metadata <- data.frame(colData(se))
pex <- assay(se)
pex[is.na(pex)] <- 0   

## PX_BINF_RvNR_W0 
PX_BINF_RvNR_W0_sample_metadata <- sample_metadata %>%
  dplyr::filter(!!sym(label_col) %in% classes,
                !!sym(timepoint_col) %in% timepoint)

PX_BINF_RvNR_W0_gex <- t(gex[,rownames(PX_BINF_RvNR_W0_sample_metadata)])

## PX_PB_RvNR_W0
PX_PB_RvNR_W0_sample_metadata <- sample_metadata %>%
  dplyr::filter(!!sym(label_col) %in% classes,
                !!sym(timepoint_col) %in% timepoint)

PX_PB_RvNR_W0_gex <- t(gex[,rownames(PX_PB_RvNR_W0_sample_metadata)])

# Save
fwrite(labels, labels_csv, quote = FALSE, row.names = FALSE)
fwrite(TX_BINF_RvNR_W0_gex, tx_binf_rvnr_w0_data_csv, quote = FALSE, row.names = FALSE)
fwrite(TX_PB_RvNR_W0_gex, tx_pb_rvnr_w0_data_csv, quote = FALSE, row.names = FALSE)
fwrite(PX_BINF_RvNR_W0_gex, px_binf_rvnr_w0_data_csv, quote = FALSE, row.names = FALSE)
fwrite(PX_PB_RvNR_W0_gex, px_pb_rvnr_w0_data_csv, quote = FALSE, row.names = FALSE)

sink(type = "message")
sink(type = "output")
close(log_file)