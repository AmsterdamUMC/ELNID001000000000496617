#!/usr/bin/env Rscript
# Supplementary Figure 1
# A) Volcanoplot OLINK responders week 8 versus baseline
# B) Volcanoplot OLINK non-responders week 8 versus baseline
# C) Scatterplot OLINK responders versus non-responders over week 8 versus baseline

# Paths
deps_R_W8vW0_csv <- snakemake@input[["degs_R_W8vW0_csv"]] #"output/260519/rna_Tissue_R_Inflamed_Week_8_vs_Baseline_deseq2_results.csv" 
deps_NR_W8vW0_csv <- snakemake@input[["degs_NR_W8vW0_csv"]] #"output/260519/rna_Tissue_NR_Inflamed_Week_8_vs_Baseline_deseq2_results.csv"

goi_xlsx <- snakemake@params[["goi_xlsx"]] #"config/features_of_interest/260716_genes.xlsx"
pwoi_xlsx <- snakemake@params[["pwoi_xlsx"]] #"config/features_of_interest/260716_genesets.xlsx"
plotting_parameters_r <- snakemake@params[["plotting_parameters"]] #"workflow/scripts/plotting_parameters.R"

figA_pdf <- snakemake@output[["figA_pdf"]]
figB_pdf <- snakemake@output[["figB_pdf"]]
figC_pdf <- snakemake@output[["figC_pdf"]]

log_file <- snakemake@log[[1]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Libraries
library(dplyr)
library(ggplot2)
library(scales)
library(ggrepel)
library(ggrastr)
library(patchwork)
library(readxl)
library(readr)
library(scales)
library(stringr)

sessionInfo()

# Import
source(plotting_parameters_r)

deps_R_W8vW0 <- read.csv(deps_R_W8vW0_csv)
deps_NR_W8vW0 <- read.csv(deps_NR_W8vW0_csv)

# Preparation

selected_features_olink_rnaseq <- c("IL6", "ICAM1", "CXCL8", "IL1B", "TNFRSF8")

## Significant results RvNR_W0

sig_deps_RvNR_W0 <- deps_RvNR_W0 %>%
  dplyr::filter(padj<0.05) %>%
  dplyr::arrange(pvalue)

## Combine degs R and NR
deps_RvNR_W8vW0_combined <- deps_R_W8vW0 %>%
  dplyr::inner_join(deps_NR_W8vW0, by = c("gene_symbol", "OlinkID", "UniProt", "Panel", "Panel_Lot_Nr"), suffix = c("_R", "_NR")) %>%
  dplyr::mutate(Significance_detailed = case_when(
    padj_R<0.05 & log2FoldChange_R<0 & pvalue_NR<0.05 & log2FoldChange_NR<0 ~ "Treatment down",
    padj_R<0.05 & log2FoldChange_R>0 & pvalue_NR<0.05 & log2FoldChange_NR>0 ~ "Treatment up",
    padj_R<0.05 & log2FoldChange_R<0 & pvalue_NR<0.05 & log2FoldChange_NR>0 ~ "R: down; NR: up",
    padj_R<0.05 & log2FoldChange_R>0 & pvalue_NR<0.05 & log2FoldChange_NR<0 ~ "R: up; NR: down",
    padj_R<0.05 & log2FoldChange_R<0 & pvalue_NR>0.05 ~ "R: down; NR: stable",
    padj_R<0.05 & log2FoldChange_R>0 & pvalue_NR>0.05 ~ "R: up; NR: stable",
    .default = "NS"
  ),
  Significance = case_when(
    padj_R<0.05 & log2FoldChange_R<0 ~ "Down",
    padj_R<0.05 & log2FoldChange_R>0 ~ "Up",
    .default = "NS"
  ),
  Significance = factor(Significance, levels = names(direction_colors)))

# A

message("Starting on A")

## Volcanoplot R:W8vW0

l2fc_R_W8vW0_range <- range(deps_RvNR_W8vW0_combined$log2FoldChange_R)

volcanoplot_deps_R_W8vW0_ggplotobj <- deps_RvNR_W8vW0_combined %>%
  dplyr::mutate(mlog10padj = -log10(padj_R),
                Significance = case_when(
                  log2FoldChange_R > 0 & padj_R < 0.05 ~ "Up",
                  log2FoldChange_R < 0 & padj_R < 0.05 ~ "Down",
                  .default = "NS")) %>%
  ggplot(aes(y = mlog10padj, x = log2FoldChange_R, color = Significance)) +
  geom_hline(yintercept = 0, color = "gray40") +
  geom_vline(xintercept = 0, color = "gray40") +
  geom_point_rast(size = 1.4, alpha = 0.8) +
  scale_color_manual(values = direction_colors,
                     drop = FALSE) +
  labs(x = "log2(fold-change)", 
       y = "-log10(padj)",
       title = "Responders: W8 vs W0") +
  guides(fill=guide_legend(title="R")) +
  # scale_y_continuous(position = "right") +
  coord_cartesian(xlim = l2fc_R_W8vW0_range) +
  theme_bw() + 
  theme(legend.position = "none",
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

## Volcanoplot NR:W8vW0

l2fc_NR_W8vW0_range <- range(deps_RvNR_W8vW0_combined$log2FoldChange_NR)

volcanoplot_deps_NR_W8vW0_ggplotobj <- deps_RvNR_W8vW0_combined %>%
  dplyr::mutate(mlog10padj = -log10(padj_NR),
                Significance = case_when(
                  log2FoldChange_NR > 0 & padj_NR < 0.05 ~ "Up",
                  log2FoldChange_NR < 0 & padj_NR < 0.05 ~ "Down",
                  .default = "NS")) %>%
  ggplot(aes(y = mlog10padj, x = log2FoldChange_NR, color = Significance)) +
  geom_hline(yintercept = 0, color = "gray40") +
  geom_vline(xintercept = 0, color = "gray40") +
  geom_point_rast(size = 1.4, alpha = 0.8) +
  scale_color_manual(values = direction_colors,
                     drop = FALSE) +
  labs(x = "log2(fold-change)", 
       y = "-log10(padj)",
       title = "Non-responders: W8 vs W0") +
  guides(fill=guide_legend(title="R")) +
  # scale_y_continuous(position = "right") +
  coord_cartesian(xlim = l2fc_NR_W8vW0_range) +
  theme_bw() + 
  theme(legend.position = "none",
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

## Scatterplot RvNR:W8vW0

cor_RvNR <- cor.test(deps_RvNR_W8vW0_combined$log2FoldChange_R, deps_RvNR_W8vW0_combined$log2FoldChange_NR)

scatterplot_deps_RvNR_W8vW0_ggplotobj <- deps_RvNR_W8vW0_combined %>%
  ggplot(aes(x = log2FoldChange_R, y = log2FoldChange_NR)) +
  geom_hline(yintercept = 0, color = "gray30") +
  geom_vline(xintercept = 0, color = "gray30") +
  geom_point_rast(data = . %>% 
                    dplyr::filter(Significance == "NS"), 
                  aes(color = Significance), size = 1.4, alpha = 0.8) +
  geom_point_rast(data = . %>% 
                    dplyr::filter(Significance != "NS"), 
                  aes(color = Significance), size = 1.4, alpha = 0.8) +
  scale_color_manual(values = direction_colors, 
                     name = "Significance") +
  scale_fill_manual(values = direction_colors,
                    drop = FALSE,
                    guide = "none") +
  annotate("text", 
           x = I(0.1), 
           y = I(0.95), 
           label = paste0("r = ", round(cor_RvNR$estimate, 3))) + 
  labs(x = "log2(W8-W0)\nResponders",
       y = "log2(W8-W0)\nNon-responders") +
  theme_bw() + 
  # scale_y_continuous(position = "right") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# Assembly

figA <- volcanoplot_deps_R_W8vW0_ggplotobj
figB <- volcanoplot_deps_NR_W8vW0_ggplotobj
figC <- scatterplot_deps_RvNR_W8vW0_ggplotobj

# Save

ggsave(filename = figA_pdf, 
       plot = figA, 
       width = 6,
       height = 6.5)
ggsave(filename = figB_pdf, 
       plot = figB, 
       width = 6,
       height = 6.5)
ggsave(filename = figC_pdf,
       plot = figC,
       width = 7,
       height = 6.5)

sink(type = "message")
sink(type = "output")
close(log_file)