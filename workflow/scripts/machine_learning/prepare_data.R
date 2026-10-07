
#!/usr/bin/env Rscript
# Export ML feature matrices + labels from the prepared RNA-seq (dds) and Olink (se) objects, for pretreatment (baseline) response prediction.

# Paths
dds_rds <- snakemake@input[["dds_rds"]]
se_rds <- snakemake@input[["se_rds"]]
mlp <- snakemake@params[["ml_params"]]
timepoint_col <- mlp$week_column
timepoint <- mlp$week
label_col <- mlp$label_column
donor_col <- mlp$donor_column
classes <- mlp$classes

labels_csv <- snakemake@output[["labels_csv"]]
tx_binf_rvnr_w0_data_csv <- snakemake@output[["tx_binf_rvnr_w0_data_csv"]]
tx_pb_rvnr_w0_data_csv <- snakemake@output[["tx_pb_rvnr_w0_data_csv"]]
px_binf_rvnr_w0_data_csv <- snakemake@output[["px_binf_rvnr_w0_data_csv"]]
px_pb_rvnr_w0_data_csv <- snakemake@output[["px_pb_rvnr_w0_data_csv"]]

log_file  <- "output/machine_learning/data/ml_prepare_data.log"

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

# Import
dds <- readRDS(dds_rds)
se <- readRDS(se_rds)
# Metadata for each modality
sample_metadata_dds <- data.frame(colData(dds))
sample_metadata_se  <- data.frame(colData(se))
# Labels
labels <- unique(sample_metadata_dds[, c(donor_col, label_col)])
labels <- labels[which(labels$Response %in% classes),]

# TX
sample_metadata <- data.frame(colData(dds))
gex <- assay(dds)

## Build the per-modality metadata subsets FIRST (before touching the assay matrices)
TX_BINF_RvNR_W0_sample_metadata <- sample_metadata_dds %>%
  dplyr::filter(.data[[label_col]] %in% classes,
                .data[[timepoint_col]] %in% timepoint,
                .data[["Tissue"]] %in% c("BINF"))
 
TX_PB_RvNR_W0_sample_metadata <- sample_metadata_dds %>%
  dplyr::filter(.data[[label_col]] %in% classes,
                .data[[timepoint_col]] %in% timepoint,
                .data[["Tissue"]] %in% c("PB"))
 
PX_BINF_RvNR_W0_sample_metadata <- sample_metadata_se %>%
  dplyr::filter(.data[[label_col]] %in% classes,
                .data[[timepoint_col]] %in% timepoint,
                .data[["Tissue"]] %in% c("BINF"))
 
PX_PB_RvNR_W0_sample_metadata <- sample_metadata_se %>%
  dplyr::filter(.data[[label_col]] %in% classes,
                .data[[timepoint_col]] %in% timepoint,
                .data[["Tissue"]] %in% c("PB"))

## Donors present in ALL four modalities
common_donors <- Reduce(
  intersect,
  list(
    TX_BINF_RvNR_W0_sample_metadata[[donor_col]],
    TX_PB_RvNR_W0_sample_metadata[[donor_col]],
    PX_BINF_RvNR_W0_sample_metadata[[donor_col]],
    PX_PB_RvNR_W0_sample_metadata[[donor_col]]
  )
)
message(sprintf("Donors shared across all 4 modalities: %d", length(common_donors)))

TX_BINF_RvNR_W0_sample_metadata <- TX_BINF_RvNR_W0_sample_metadata[TX_BINF_RvNR_W0_sample_metadata[[donor_col]] %in% common_donors, , drop = FALSE]
TX_PB_RvNR_W0_sample_metadata  <- TX_PB_RvNR_W0_sample_metadata[TX_PB_RvNR_W0_sample_metadata[[donor_col]] %in% common_donors, , drop = FALSE]
PX_BINF_RvNR_W0_sample_metadata <- PX_BINF_RvNR_W0_sample_metadata[PX_BINF_RvNR_W0_sample_metadata[[donor_col]] %in% common_donors, , drop = FALSE]
PX_PB_RvNR_W0_sample_metadata  <- PX_PB_RvNR_W0_sample_metadata[PX_PB_RvNR_W0_sample_metadata[[donor_col]] %in% common_donors, , drop = FALSE]
labels <- labels[labels[[donor_col]] %in% common_donors, , drop = FALSE]
# TX
gex <- assay(dds)
 
#vsd <- DESeq2::vst(dds)
#gex <- assay(vsd)
 
## TX_BINF_RvNR_W0
TX_BINF_RvNR_W0_gex <- t(gex[, rownames(TX_BINF_RvNR_W0_sample_metadata), drop = FALSE])
TX_BINF_RvNR_W0_keep <- colSums(TX_BINF_RvNR_W0_gex >= 10, na.rm = TRUE) >= 5
TX_BINF_RvNR_W0_gex <- TX_BINF_RvNR_W0_gex[, TX_BINF_RvNR_W0_keep, drop = FALSE]
TX_BINF_RvNR_W0_depth <- rowSums(TX_BINF_RvNR_W0_gex, na.rm = TRUE)
TX_BINF_RvNR_W0_gex <- TX_BINF_RvNR_W0_gex / (rowSums(TX_BINF_RvNR_W0_gex) / 1e6)
TX_BINF_RvNR_W0_gex <- log2(TX_BINF_RvNR_W0_gex + 1)
 
## TX_PB_RvNR_W0
TX_PB_RvNR_W0_gex <- t(gex[, rownames(TX_PB_RvNR_W0_sample_metadata), drop = FALSE])
TX_PB_RvNR_W0_keep <- colSums(TX_PB_RvNR_W0_gex >= 10, na.rm = TRUE) >= 5
TX_PB_RvNR_W0_gex <- TX_PB_RvNR_W0_gex[, TX_PB_RvNR_W0_keep, drop = FALSE]
TX_PB_RvNR_W0_depth <- rowSums(TX_PB_RvNR_W0_gex, na.rm = TRUE)
TX_PB_RvNR_W0_gex <- TX_PB_RvNR_W0_gex / (rowSums(TX_PB_RvNR_W0_gex) / 1e6)
TX_PB_RvNR_W0_gex <- log2(TX_PB_RvNR_W0_gex + 1)
 
# Olink
pex <- assay(se)
pex[is.na(pex)] <- 0
 
## PX_BINF_RvNR_W0
PX_BINF_RvNR_W0_pex <- t(pex[, rownames(PX_BINF_RvNR_W0_sample_metadata), drop = FALSE])
 
## PX_PB_RvNR_W0
PX_PB_RvNR_W0_pex <- t(pex[, rownames(PX_PB_RvNR_W0_sample_metadata), drop = FALSE])

rownames(TX_BINF_RvNR_W0_gex)  <- TX_BINF_RvNR_W0_sample_metadata[[donor_col]]
rownames(TX_PB_RvNR_W0_gex)    <- TX_PB_RvNR_W0_sample_metadata[[donor_col]]
rownames(PX_BINF_RvNR_W0_pex)  <- PX_BINF_RvNR_W0_sample_metadata[[donor_col]]
rownames(PX_PB_RvNR_W0_pex)    <- PX_PB_RvNR_W0_sample_metadata[[donor_col]]

# Save
fwrite(labels, labels_csv, quote = FALSE, row.names = FALSE)
fwrite(TX_BINF_RvNR_W0_gex, tx_binf_rvnr_w0_data_csv, quote = FALSE, row.names = TRUE)
fwrite(TX_PB_RvNR_W0_gex, tx_pb_rvnr_w0_data_csv, quote = FALSE, row.names = TRUE)
fwrite(PX_BINF_RvNR_W0_pex, px_binf_rvnr_w0_data_csv, quote = FALSE, row.names = TRUE)
fwrite(PX_PB_RvNR_W0_pex, px_pb_rvnr_w0_data_csv, quote = FALSE, row.names = TRUE)
 
sink(type = "message")
sink(type = "output")
close(log_file)