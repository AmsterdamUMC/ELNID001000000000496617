#!/usr/bin/env Rscript
# Download MSigDbR Hallmark, KEGG, and Reactome gene sets

# Libraries
library(msigdbr)
library(tidyverse)

# Paths
gs_msigdb_rds <- snakemake@output[["gs_rds"]]
log_file <- snakemake@log[[1]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Load Hallmark collection
message("Loading Hallmark collection")
msig_hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::distinct(gs_name, gene_symbol)

# Load C2 collection with REACTOME
message("Loading C2 collection REACTOME")
msig_reactome <- msigdbr(species = "Homo sapiens", collection = "C2") %>%
  dplyr::filter(gs_subcollection %in% c("CP:REACTOME")) %>%
  dplyr::distinct(gs_name, gene_symbol)

# Load C2 collection with KEGG
message("Loading C2 collection KEGG")
msig_kegg <- msigdbr(species = "Homo sapiens", collection = "C2") %>%
  dplyr::filter(gs_subcollection %in% c("CP:KEGG_LEGACY")) %>%
  dplyr::distinct(gs_name, gene_symbol)

# Combine 
msig <- dplyr::bind_rows(msig_hallmark, msig_reactome, msig_kegg)
gs_msigdb <- split(msig$gene_symbol, msig$gs_name)

# Save
saveRDS(gs_msigdb, gs_msigdb_rds, compress = "gzip")

sessionInfo()

sink(type = "message")
sink(type = "output")
close(log_file)
