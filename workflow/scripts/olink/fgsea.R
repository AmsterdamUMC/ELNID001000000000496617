#!/usr/bin/env Rscript
# Perform the gene set enrichment analysis per comparison

# Paths
deps_csv <- snakemake@input[["deps_csv"]]
msigdb_rds <- snakemake@input[["msigdb_rds"]]
fgsea_params <- snakemake@params[["fgsea_params"]]
log_file <- snakemake@log[[1]]
fgsea_csv <- snakemake@output[["fgsea_csv"]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Libraries
library(tidyverse)
library(fgsea)
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)

sessionInfo()

# Import full dds object
deps <- read.csv(deps_csv)

gs_msigdb <- readRDS(msigdb_rds)

deps_rank <- deps %>%
  filter(!is.na(gene_symbol)) %>%
  group_by(gene_symbol) %>%
  slice_max(order_by = stat, n = 1) %>%
  summarise(ranks = stat[1]) %>%
  ungroup()
ranks <- deps_rank$ranks
names(ranks) <- deps_rank$gene_symbol
  
ranks <- sort(ranks, decreasing = TRUE)
print(fgsea_params)
fgsea_results <- fgseaMultilevel(pathways = gs_msigdb,
                                 stats = ranks,
                                 minSize = fgsea_params$min_size,
                                 maxSize = fgsea_params$max_size,
                                 nPermSimple = fgsea_params$n_perm_simple)
fgsea_results <- fgsea_results[order(fgsea_results$pval),]

# Save
fwrite(fgsea_results, fgsea_csv)

sink(type = "message")
sink(type = "output")
close(log_file)