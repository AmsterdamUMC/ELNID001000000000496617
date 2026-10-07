#!/usr/bin/env Rscript
# Supplementary Figure 1
# A) Scatterplot of the differential gene expression 
# B) Scatterplot of the differential gene sets 
# C) Arrowplot of selected pathways

# Paths
degs_NR_W8vW0_TOF_csv <- snakemake@input[["degs_NR_W8vW0_TOF_csv"]]
fgsea_R_W8vW0_TOF_csv <- snakemake@input[["fgsea_R_W8vW0_TOF_csv"]]
fgsea_NR_W8vW0_TOF_csv <- snakemake@input[["fgsea_NR_W8vW0_TOF_csv"]]
degs_NR_W8vW0_TOF_melonardanaz2025_csv <- snakemake@input[["degs_NR_W8vW0_TOF_melonardanaz2025_csv"]]
fgsea_R_W8vW0_TOF_melonardanaz2025_csv <- snakemake@input[["fgsea_R_W8vW0_TOF_melonardanaz2025_csv"]]
fgsea_NR_W8vW0_TOF_melonardanaz2025_csv <- snakemake@input[["fgsea_NR_W8vW0_TOF_melonardanaz2025_csv"]]

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

degs_NR_W8vW0_TOF <- read.csv(degs_NR_W8vW0_TOF_csv)
fgsea_R_W8vW0_TOF <- read.csv(fgsea_R_W8vW0_TOF_csv)
fgsea_NR_W8vW0_TOF <- read.csv(fgsea_NR_W8vW0_TOF_csv)
degs_NR_W8vW0_TOF_melonardanaz2025 <- read.csv(degs_NR_W8vW0_TOF_melonardanaz2025_csv)
fgsea_R_W8vW0_TOF_melonardanaz2025 <- read.csv(fgsea_R_W8vW0_TOF_melonardanaz2025_csv)
fgsea_NR_W8vW0_TOF_melonardanaz2025 <- read.csv(fgsea_NR_W8vW0_TOF_melonardanaz2025_csv)

# Preparation

pwoi_figure4 <- pwoi %>%
  dplyr::filter(!is.na(Figure4))

## Combine degs R and NR
degs_NR_W8vW0_TOF_ownvmelonardanaz2025 <- degs_NR_W8vW0_TOF %>%
  dplyr::inner_join(degs_NR_W8vW0_TOF_melonardanaz2025, by = c("gene_symbol", "ENSG"), suffix = c("_own", "_melonardanaz2025")) %>%
  dplyr::mutate(Significance_detailed = case_when(
    padj_own<0.05 & log2FoldChange_own<0 & pvalue_melonardanaz2025<0.05 & log2FoldChange_melonardanaz2025<0 ~ "Consistent down",
    padj_own<0.05 & log2FoldChange_own>0 & pvalue_melonardanaz2025<0.05 & log2FoldChange_melonardanaz2025>0 ~ "Consistent up",
    padj_own<0.05 & log2FoldChange_own<0 & pvalue_melonardanaz2025<0.05 & log2FoldChange_melonardanaz2025>0 ~ "Inconsistent",
    padj_own<0.05 & log2FoldChange_own>0 & pvalue_melonardanaz2025<0.05 & log2FoldChange_melonardanaz2025<0 ~ "Inconsistent",
    padj_own<0.05 & log2FoldChange_own<0 & pvalue_melonardanaz2025>0.05 ~ "Inconsistent",
    padj_own<0.05 & log2FoldChange_own>0 & pvalue_melonardanaz2025>0.05 ~ "Inconsistent",
    .default = "NS"
  ),
  Significance = case_when(
    Significance_detailed == "Consistent down" ~ "Down",
    Significance_detailed == "Consistent up" ~ "Up",
    .default = "NS"
  ),
  Significance = factor(Significance, levels = names(direction_colors)))

## Combine fgsea NR wide
fgsea_NR_W8vW0_TOF_ownvmelonardanaz2025_wide <- fgsea_NR_W8vW0_TOF %>%
  dplyr::select(pathway, pval, padj, NES) %>%
  dplyr::inner_join(fgsea_NR_W8vW0_TOF_melonardanaz2025 %>%
                       dplyr::select(pathway, pval, padj, NES),
                    by = "pathway", suffix = c("_own", "_melonardanaz2025")) %>%
  dplyr::mutate(label = gsub("_", " ", pathway),
                label = gsub("(REACTOME|HALLMARK|KEGG)(.+)", "(\\1)\\2", label)) %>%
  dplyr::mutate(Significance = case_when(
    padj_own<0.05 & NES_own<0 & padj_melonardanaz2025<0.05 & NES_melonardanaz2025<0 ~ "Consistent down",
    padj_own<0.05 & NES_own>0 & padj_melonardanaz2025<0.05 & NES_melonardanaz2025>0 ~ "Consistent up",
    padj_own<0.05 & NES_own<0 & padj_melonardanaz2025<0.05 & NES_melonardanaz2025>0 ~ "Own: down; MA2025: up",
    padj_own<0.05 & NES_own>0 & padj_melonardanaz2025<0.05 & NES_melonardanaz2025<0 ~ "Own: up; MA2025: down",
    padj_own<0.05 & NES_own<0 & padj_melonardanaz2025>0.05 ~ "Own: up",
    padj_own<0.05 & NES_own>0 & padj_melonardanaz2025>0.05 ~ "Own: down",
    padj_own>0.05 & padj_melonardanaz2025<0.05 & NES_melonardanaz2025>0 ~ "MA2025: up",
    padj_own>0.05 & padj_melonardanaz2025<0.05 & NES_melonardanaz2025<0 ~ "MA2025: down",
    .default = "NS"
  ),
  Significance = factor(Significance, levels = rev(c("Consistent up", "Consistent down", "Own: up", "Own: down", "MA2025: up", "MA2025: down", "Own: up; MA2025: down", "Own: down; MA2025: up", "NS"))))

## Combine fgsea R and NR long
fgsea_RvNR_W8vW0_TOF_ownvmelonardanaz2025_long <- fgsea_R_W8vW0_TOF %>%
  dplyr::select(pathway, pval, padj, NES) %>%
  dplyr::mutate(Treatment = "Tofacitinib",
                Source = "Own",
                Response = "R") %>%
  dplyr::rows_append(fgsea_NR_W8vW0_TOF %>%
                       dplyr::select(pathway, pval, padj, NES) %>%
                       dplyr::mutate(Treatment = "Tofacitinib",
                                     Source = "Own",
                                     Response = "NR")) %>%
  dplyr::rows_append(fgsea_R_W8vW0_TOF_melonardanaz2025 %>%
                       dplyr::select(pathway, pval, padj, NES) %>%
                       dplyr::mutate(Treatment = "Tofacitinib",
                                     Source = "MA2025",
                                     Response = "R")) %>%
  dplyr::rows_append(fgsea_NR_W8vW0_TOF_melonardanaz2025 %>%
                       dplyr::select(pathway, pval, padj, NES) %>%
                       dplyr::mutate(Treatment = "Tofacitinib",
                                     Source = "MA2025",
                                     Response = "NR")) %>%
  dplyr::mutate(label = gsub("_", " ", pathway),
                label = gsub("(REACTOME|HALLMARK|KEGG)(.+)", "(\\1)\\2", label)) %>%
  dplyr::mutate(Direction = ifelse(NES<0, "Down", "Up"),
                Significance = ifelse(padj<0.05, "Significant", "NS"),
                Response = factor(Response, levels = c("NR", "R")),
                Cohort = factor(paste0(Treatment, "\n", Source), levels = c("Tofacitinib\nOwn", "Tofacitinib\nMA2025", "Infliximab\nT2011")))
# A

message("Starting on A")

## Scatterplot NR:W8vW0 Own vs Melon-Ardanaz 2025

cor_degs_NR_W8vW0_TOF_ownvmelonardanaz2025 <- cor.test(degs_NR_W8vW0_TOF_ownvmelonardanaz2025$log2FoldChange_own, degs_NR_W8vW0_TOF_ownvmelonardanaz2025$log2FoldChange_melonardanaz2025)

scatterplot_degs_NR_W8vW0_ggplotobj <- degs_NR_W8vW0_TOF_ownvmelonardanaz2025 %>%
  ggplot(aes(x = log2FoldChange_own, y = log2FoldChange_melonardanaz2025)) +
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
           label = paste0("r = ", round(cor_degs_NR_W8vW0_TOF_ownvmelonardanaz2025$estimate, 3))) + 
  labs(title = "Tofacitinib",
       subtitle = "NR: W8 - W0",
       x = "log2(w8-w0)\nOwn",
       y = "log2(w8-w0)\nMelon-Ardanaz 2025") +
  theme_bw() + 
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# B

cor_fgsea_NR_W8vW0_TOF_ownvmelonardanaz2025 <- cor.test(fgsea_NR_W8vW0_TOF_ownvmelonardanaz2025_wide$NES_own, fgsea_NR_W8vW0_TOF_ownvmelonardanaz2025_wide$NES_melonardanaz2025)

scatterplot_fgsea_NR_W8vW0_ggplotobj <- fgsea_NR_W8vW0_TOF_ownvmelonardanaz2025_wide %>%
  ggplot(aes(x = NES_own, y = NES_melonardanaz2025)) +
  geom_hline(yintercept = 0, color = "gray30") +
  geom_vline(xintercept = 0, color = "gray30") +
  geom_point_rast(data = . %>% 
                    dplyr::filter(Significance == "NS"), 
                  aes(color = Significance), size = 1.4, alpha = 0.8) +
  geom_point_rast(data = . %>% 
                    dplyr::filter(Significance != "NS"), 
                  aes(color = Significance), size = 1.4, alpha = 0.8) +
  scale_color_manual(values = ownvma2025_direction_colors,
  name = "Significance") +
  scale_fill_manual(values = ownvma2025_direction_colors,
  drop = FALSE,
  guide = "none") +
  annotate("text", 
           x = I(0.1), 
           y = I(0.95), 
           label = paste0("r = ", round(cor_fgsea_NR_W8vW0_TOF_ownvmelonardanaz2025$estimate, 3))) + 
  labs(title = "Tofacitinib",
       subtitle = "NR: W8 - W0",
       x = "NES\nOwn",
       y = "NES\nMelon-Ardanaz 2025") +
  theme_bw() + 
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# C

arrowplot_fgsea_RvNR_W8vW0_TOF_ownvmelonardanaz2025_ggplotobj <- fgsea_RvNR_W8vW0_TOF_ownvmelonardanaz2025_long %>%
  dplyr::inner_join(pwoi_figure4, by = c("pathway"= "Geneset")) %>%
  ggplot(aes(y = Response, x = Label)) +
  geom_point(aes(fill = NES, size = -log10(padj), alpha = Significance, shape = Direction)) +
  scale_shape_manual(values = direction_shapes) +
  scale_alpha_manual(values = significance_alpha) +
  scale_fill_gradient2(high = "coral3", low = "deepskyblue3", mid = "white") +
  facet_grid(Cohort~Figure4, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

# Assembly

figA <- scatterplot_degs_NR_W8vW0_ggplotobj
figB <- scatterplot_fgsea_NR_W8vW0_ggplotobj
figC <- arrowplot_fgsea_RvNR_W8vW0_TOF_ownvmelonardanaz2025_ggplotobj

# Save

ggsave(filename = figA_pdf, 
       plot = figA, 
       width = 6,
       height = 6.5)
ggsave(filename = figB_pdf, 
       plot = figB, 
       width = 7.5,
       height = 6.5)
ggsave(filename = figC_pdf,
       plot = figC,
       width = 8,
       height = 6)

sink(type = "message")
sink(type = "output")
close(log_file)