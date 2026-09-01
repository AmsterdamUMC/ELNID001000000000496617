#!/usr/bin/env R
# Scatterplot of the differentially expressed genes effect sizes (log2foldchange) obtained from comparing the week 8 with week 0 from  responders comparing todacitinib with infliximab (Toedtger) (DOI: XXX)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste0("Script needs 5 arguments. Current input is:", args))
}

# Libraries
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrastr)
library(readxl)
library(readr)

# Paths

## In
degs_BINF_R_w8vw0_TOF_csv <- args[1] #"output/260519/Infliximab_R_Week_8_vs_Baseline_deseq2_results.csv"
degs_BINF_NR_w8vw0_TOF_csv <- args[2] #"output/260519/Infliximab_NR_Week_8_vs_Baseline_deseq2_results.csv"

## Out
scatterplot_degs_BINF_RvNR_w8vw0_combined_pdf <- args[5]

# Import
degs_BINF_R_w8vw0_TOF <- read.csv(degs_BINF_R_w8vw0_TOF_csv)
degs_BINF_NR_w8vw0_TOF <- read.csv(degs_BINF_NR_w8vw0_TOF_csv)
degs_BINF_R_w8vw0_IFX <- read_excel(degs_BINF_R_w8vw0_IFX_xlsx)
degs_BINF_NR_w8vw0_IFX <- read_excel(degs_BINF_NR_w8vw0_IFX_xlsx)

# Combine results R and NR
degs_BINF_R_w8vw0_combined <- degs_BINF_R_w8vw0_TOF %>%
  dplyr::inner_join(degs_BINF_NR_w8vw0_TOF, by = c("gene_symbol" = "gene"), suffix = c("_own", "_ma2025")) %>%
  dplyr::mutate(Response = "R")
degs_BINF_NR_w8vw0_combined <- degs_BINF_NR_w8vw0_TOF %>%
  dplyr::inner_join(degs_BINF_NR_w8vw0_IFX, by = c("gene_symbol" = "gene"), suffix = c("_own", "_ma2025")) %>%
  dplyr::mutate(Response = "NR")

degs_BINF_RvNR_w8vw0_combined <- degs_BINF_R_w8vw0_combined %>%
  dplyr::rows_append(degs_BINF_NR_w8vw0_combined) %>%
  dplyr::mutate(Response = factor(Response, levels = c("R", "NR")))

# Plot
scatterplot_degs_BINF_RvNR_w8vw0_combined_ggplotobj <- degs_BINF_RvNR_w8vw0_combined %>%
  ggplot(aes(x = log2FoldChange_own, y = log2FoldChange_ma2025)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point_rast() +
  geom_smooth(method = "lm") +
  labs(title = "Comparison with Ardanaz Melon 2025\nDOI: 10.1093/ecco-jcc/jjaf076",
       subtitle = "log2 (fold change)",
       x = "Own\n(bulk RNAseq)",
       y = "Ardanaz Melon 2025\n(pseudobulk scRNAseq)") +
  facet_grid(.~Response) +
  theme_bw() +
  xlim(c(-10, 10)) +
  ylim(c(-10, 10)) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())

# Save
ggsave(filename = scatterplot_degs_BINF_RvNR_w8vw0_combined_pdf, 
       plot = scatterplot_degs_BINF_RvNR_w8vw0_combined_ggplotobj, 
       width = 9,
       height = 5)