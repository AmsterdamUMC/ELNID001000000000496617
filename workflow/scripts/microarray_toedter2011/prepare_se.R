#!/usr/bin/env Rscript
# Build the general SummarizedExperiment from the microarray expression values and sample metadata as stored on GSE23597

# Paths
toedter2011_gse <- snakemake@input[["toedter2011_gse"]]
se_rds <- snakemake@output[["se_rds"]]
log_file <- snakemake@log[[1]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Libraries
library(tidyverse)
library(GEOquery)
library(SummarizedExperiment)

sessionInfo()

gse <- getGEO(filename = toedter2011_gse)

# Prepare sample metadata

sample_metadata <- pData(gse) %>%
  dplyr::rename(SampleID = geo_accession,
                DonorID = `characteristics_ch1.3`,
                Dosage = `dose:ch1`,
                Response = `wk8 response:ch1`,
                Week = `time:ch1`) %>%
  dplyr::select(SampleID, DonorID, Dosage, Response, Week) %>%
  dplyr::mutate(DonorID = gsub("subject: ", "", DonorID),
                Treatment = ifelse(Dosage %in% c("5mg/kg", "10mg/kg"), "Infliximab", "Placebo"),
                Dosage = gsub("([0-9]+)mg/kg", "\\1", Dosage),
                Dosage = as.numeric(ifelse(Dosage == "placebo", NA, Dosage)),
                Response = factor(ifelse(Response == "No", "NR", "R"), levels = c("NR", "R")),
                Tissue = "BINF",
                Week = factor(Week, levels = c("W0", "W8", "W30")))

message("--- Creating SummarizedExperiment ---")
toedter_se <- SummarizedExperiment(assays = list(exprs = exprs(gse),
                                                 log = log2(exprs(gse)+1)),
                                   rowData = fData(gse),
                                   colData = sample_metadata)

# Save

saveRDS(toedter_se, se_rds, compress = "gzip")
message("Saved: ", se_rds)

sink(type = "output")
sink(type = "message")
close(log_file)