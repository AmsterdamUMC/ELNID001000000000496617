#!/usr/bin/env Rscript
# Figure 4
# A) Correlation overlap between bulk RNA-sequencing and pseudobulked single-cell RNA-sequencing.
# B) Correlation overlap between tofacitinib and infliximab.
# C) 
# D) Volcanoplots representing the longitudinal differentially expressed genes identified when comparing week 8 with week 0 in responders (top-right) and non-responders (bottom-left) as well as a scatterplot (bottom-right) representing the longitudinal log2fold changes for responders (X-axis) and non-responders (Y-axis).
# E)

# Paths
vst_R_W8vW0_TOF_rds <- snakemake@input[["vst_R_W8vW0_TOF_rds"]] #"output/rnaseq/subsets/dds_general.rds"
microarray_toedter2011_R_W8vW0_IFX_se_rds <- snakemake@input[["microarray_toedter2011_R_W8vW0_IFX_se_rds"]] # "output/microarray/GSE23597_series_matrix.txt"
degs_R_W8vW0_TOF_csv <- snakemake@input[["degs_R_W8vW0_TOF_csv"]] #"output/260519/rna_Tissue_R_Inflamed_Week_8_vs_Baseline_deseq2_results.csv"
fgsea_R_W8vW0_TOF_csv <- snakemake@input[["fgsea_R_W8vW0_TOF_csv"]] #"output/260519/rna_Tissue_R_Inflamed_Week_8_vs_Baseline_fgsea_results.csv"
fgsea_NR_W8vW0_TOF_csv <- snakemake@input[["fgsea_NR_W8vW0_TOF_csv"]] #"output/260519/rna_Tissue_NR_Inflamed_Week_8_vs_Baseline_fgsea_results.csv"
degs_R_W8vW0_TOF_melonardanaz2025_csv <- snakemake@input[["degs_R_W8vW0_TOF_melonardanaz2025_csv"]] #"output/260615/rna_Tissue_R_Inflamed_Week_8_vs_Baseline_deseq2_results_melonardanaz2025.csv"
fgsea_R_W8vW0_TOF_melonardanaz2025_csv <- snakemake@input[["fgsea_R_W8vW0_TOF_melonardanaz2025_csv"]] #"output/260615/rna_Tissue_R_Inflamed_Week_8_vs_Baseline_fgsea_results_melonardanaz2025.csv"
fgsea_NR_W8vW0_TOF_melonardanaz2025_csv <- snakemake@input[["fgsea_NR_W8vW0_TOF_melonardanaz2025_csv"]] #"output/260615/rna_Tissue_NR_Inflamed_Week_8_vs_Baseline_fgsea_results_melonardanaz2025.csv"
degs_R_W8vW0_IFX_toedter2011_csv <- snakemake@input[["degs_R_W8vW0_IFX_toedter2011_csv"]] #"output/260519/Infliximab_R_Week_8_vs_Baseline_limma_results.csv"
fgsea_R_W8vW0_IFX_toedter2011_csv <- snakemake@input[["fgsea_R_W8vW0_IFX_toedter2011_csv"]] #"output/260519/Infliximab_R_Week_8_vs_Baseline_fgsea_results.csv"
fgsea_NR_W8vW0_IFX_toedter2011_csv <- snakemake@input[["fgsea_NR_W8vW0_IFX_toedter2011_csv"]] #"output/260519/Infliximab_NR_Week_8_vs_Baseline_fgsea_results.csv"

goi_xlsx <- snakemake@params[["goi_xlsx"]] #"config/features_of_interest/260716_genes.xlsx"
pwoi_xlsx <- snakemake@params[["pwoi_xlsx"]] #"config/features_of_interest/260716_genesets.xlsx"
plotting_parameters_r <- snakemake@params[["plotting_parameters"]] #"workflow/scripts/plotting_parameters.R"

figA_pdf <- snakemake@output[["figA_pdf"]]
figBD_pdf <- snakemake@output[["figBD_pdf"]]
figC_pdf <- snakemake@output[["figC_pdf"]]
figE_pdf <- snakemake@output[["figE_pdf"]]

log_file <- snakemake@log[[1]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Libraries
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrastr)
library(patchwork)
library(readxl)
library(readr)
library(stringr)
library(ComplexHeatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Import
source(plotting_parameters_r)
goi <- readxl::read_excel(goi_xlsx)
pwoi <- readxl::read_excel(pwoi_xlsx)

vst_R_W8vW0_TOF <- readRDS(vst_R_W8vW0_TOF_rds)
print(degs_R_W8vW0_TOF_csv)
degs_R_W8vW0_TOF <- read.csv(degs_R_W8vW0_TOF_csv)
print(fgsea_R_W8vW0_TOF_csv)
fgsea_R_W8vW0_TOF <- read.csv(fgsea_R_W8vW0_TOF_csv)
print(fgsea_NR_W8vW0_TOF_csv)
fgsea_NR_W8vW0_TOF <- read.csv(fgsea_NR_W8vW0_TOF_csv)
print(degs_R_W8vW0_TOF_melonardanaz2025_csv)
degs_R_W8vW0_TOF_melonardanaz2025 <- read.csv(degs_R_W8vW0_TOF_melonardanaz2025_csv)
print(fgsea_R_W8vW0_TOF_melonardanaz2025_csv)
fgsea_R_W8vW0_TOF_melonardanaz2025 <- read.csv(fgsea_R_W8vW0_TOF_melonardanaz2025_csv)
print(fgsea_NR_W8vW0_TOF_melonardanaz2025_csv)
fgsea_NR_W8vW0_TOF_melonardanaz2025 <- read.csv(fgsea_NR_W8vW0_TOF_melonardanaz2025_csv)

microarray_toedter2011_R_W8vW0_IFX_se <- readRDS(microarray_toedter2011_R_W8vW0_IFX_se_rds)
print(degs_R_W8vW0_IFX_toedter2011_csv)
degs_R_W8vW0_IFX_toedter2011 <- read.csv(degs_R_W8vW0_IFX_toedter2011_csv)
print(fgsea_R_W8vW0_IFX_toedter2011_csv)
fgsea_R_W8vW0_IFX_toedter2011 <- read.csv(fgsea_R_W8vW0_IFX_toedter2011_csv)
print(fgsea_NR_W8vW0_IFX_toedter2011_csv)
fgsea_NR_W8vW0_IFX_toedter2011 <- read.csv(fgsea_NR_W8vW0_IFX_toedter2011_csv)

# Preparation

goi <- goi %>%
  dplyr::mutate(ENSG = mapIds(org.Hs.eg.db, 
                              keys = Gene,
                              column = "ENSEMBL", 
                              keytype = "SYMBOL"))

goi_figure4 <- goi %>%
  dplyr::filter(!is.na(Figure4))

pwoi_figure4 <- pwoi %>%
  dplyr::filter(!is.na(Figure4))

## Combine degs own and Melon-Ardanaz 2025

degs_R_W8vW0_TOF_ownvmelonardanaz2025 <- degs_R_W8vW0_TOF %>%
  dplyr::inner_join(degs_R_W8vW0_TOF_melonardanaz2025, by = c("gene_symbol", "ENSG"), suffix = c("_own", "_melonardanaz2025")) %>%
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

# Combine degs tofacitinib and infliximab
degs_R_W8vW0_TOFvIFX <- degs_R_W8vW0_TOF %>%
  dplyr::inner_join(degs_R_W8vW0_IFX_toedter2011 %>%
                      dplyr::group_by(gene_symbol) %>%
                      dplyr::summarize(log2FoldChange = mean(log2FoldChange),
                                       AveExpr = mean(AveExpr),
                                       stat = mean(stat),
                                       pvalue = mean(pvalue)) %>%
                      dplyr::mutate(padj = p.adjust(pvalue, method = "BH")), 
                    by = c("gene_symbol"), 
                    suffix = c("_TOF", "_IFX")) %>%
  dplyr::mutate(Significance_detailed = case_when(
    padj_TOF<0.05 & log2FoldChange_TOF<0 & pvalue_IFX<0.05 & log2FoldChange_IFX<0 ~ "Treatment down",
    padj_TOF<0.05 & log2FoldChange_TOF>0 & pvalue_IFX<0.05 & log2FoldChange_IFX>0 ~ "Treatment up",
    padj_TOF<0.05 & log2FoldChange_TOF<0 & pvalue_IFX<0.05 & log2FoldChange_IFX>0 ~ "TOF: down; IFX: up",
    padj_TOF<0.05 & log2FoldChange_TOF>0 & pvalue_IFX<0.05 & log2FoldChange_IFX<0 ~ "TOF: up; IFX: down",
    padj_TOF<0.05 & log2FoldChange_TOF<0 & pvalue_IFX>0.05 ~ "TOF: down; IFX: stable",
    padj_TOF<0.05 & log2FoldChange_TOF>0 & pvalue_IFX>0.05 ~ "TOF: up; IFX: stable",
    .default = "NS"
  ))

degs_R_W8vW0_TOFvIFX %>% 
  dplyr::left_join(goi_figure4, by = c("gene_symbol" = "Gene")) %>%
  dplyr::filter(!is.na(Category)) %>%
  dplyr::select(gene_symbol, Category, log2FoldChange_TOF, log2FoldChange_IFX, pvalue_TOF, pvalue_IFX, Significance_detailed) %>%
  dplyr::arrange(Category, Significance_detailed)

degs_R_W8vW0_TOFvIFX %>% 
  dplyr::left_join(goi, by = c("gene_symbol" = "Gene")) %>%
  dplyr::filter(!is.na(Category)) %>%
  dplyr::select(gene_symbol, Category, log2FoldChange_TOF, log2FoldChange_IFX, pvalue_TOF, pvalue_IFX, Significance_detailed) %>%
  dplyr::arrange(Category, Significance_detailed) %>%
  dplyr::filter(Significance_detailed == "Treatment down")

degs_R_W8vW0_TOFvIFX %>% 
  dplyr::left_join(goi, by = c("gene_symbol" = "Gene")) %>%
  dplyr::filter(!is.na(Category)) %>%
  dplyr::select(gene_symbol, Category, log2FoldChange_TOF, log2FoldChange_IFX, pvalue_TOF, pvalue_IFX, Significance_detailed) %>%
  dplyr::arrange(Category, Significance_detailed) %>%
  dplyr::filter(!Significance_detailed %in% c("Treatment down", "NS"))

# Combine fgsea tofacitinib own, tofacitinib Melon-Ardanaz 2025 and infliximab Toedter 2011
fgsea_RvNR_W8vW0_TOFvIFX_ownvmelonardanaz2025vtoedter2011 <- fgsea_R_W8vW0_TOF %>%
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
  dplyr::rows_append(fgsea_R_W8vW0_IFX_toedter2011 %>%
                       dplyr::select(pathway, pval, padj, NES) %>%
                       dplyr::mutate(Treatment = "Infliximab",
                                     Source = "T2011",
                                     Response = "R")) %>%
  dplyr::rows_append(fgsea_NR_W8vW0_IFX_toedter2011 %>%
                       dplyr::select(pathway, pval, padj, NES) %>%
                       dplyr::mutate(Treatment = "Infliximab",
                                     Source = "T2011",
                                     Response = "NR")) %>%
  dplyr::mutate(label = gsub("_", " ", pathway),
                label = gsub("(REACTOME|HALLMARK|KEGG)(.+)", "(\\1)\\2", label)) %>%
  dplyr::mutate(Direction = ifelse(NES<0, "Down", "Up"),
                Significance = ifelse(padj<0.05, "Significant", "NS"),
                Response = factor(Response, levels = c("NR", "R")),
                Cohort = factor(paste0(Treatment, "\n", Source), levels = c("Tofacitinib\nOwn", "Tofacitinib\nMA2025", "Infliximab\nT2011")))
  
# A

cor_ownvmelonardanaz2025 <- cor.test(degs_R_W8vW0_TOF_ownvmelonardanaz2025$log2FoldChange_own, degs_R_W8vW0_TOF_ownvmelonardanaz2025$log2FoldChange_melonardanaz2025)

scatterplot_degs_r_W8vW0_ownvmelonardanaz2025_ggplotobj <- degs_R_W8vW0_TOF_ownvmelonardanaz2025 %>%
  ggplot(aes(x = log2FoldChange_own, y = log2FoldChange_melonardanaz2025)) +
  geom_hline(yintercept = 0, color = "gray30") +
  geom_vline(xintercept = 0, color = "gray30") +
  geom_point_rast(data = . %>%
                    dplyr::filter(Significance == "NS"),
                  aes(color = Significance)) +
  geom_point_rast(data = . %>%
                    dplyr::filter(Significance != "NS"),
                  aes(color = Significance)) +
  scale_color_manual(values = direction_colors) +
  annotate("text", 
           x = I(0.1), 
           y = I(0.95), 
           label = paste0("r = ", round(cor_ownvmelonardanaz2025$estimate, 3))) + 
  labs(title = "Tofacitinib validation",
       subtitle = "W8 - W0",
       x = "log2(w8-w0)\nOwn",
       y = "log2(w8-w0)\nMelon-Ardanaz 2025") +
  theme_bw() + 
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# B/D v1

arrowplot_fgsea_RvNR_W8vW0_TOFvIFX_ownvmelonardanaz2025vtoedter2011_ggplotobj <- fgsea_RvNR_W8vW0_TOFvIFX_ownvmelonardanaz2025vtoedter2011 %>%
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

# B/D v2

arrowplot_fgsea_RvNR_W8vW0_TOFvIFX_ownvmelonardanaz2025vtoedter2011_ggplotobj <- fgsea_RvNR_W8vW0_TOFvIFX_ownvmelonardanaz2025vtoedter2011 %>%
  dplyr::inner_join(pwoi_figure4, by = c("pathway"= "Geneset")) %>%
  dplyr::filter(Response == "R") %>%
  dplyr::mutate(Cohort = factor(Cohort, levels = rev(levels(Cohort))),
                Drug = factor(ifelse(Treatment == "Tofacitinib", "TOF", "IFX"), levels = c("TOF", "IFX")),
                Source = case_when(
                  Cohort == "Tofacitinib\nOwn" ~ "Own",
                  Cohort == "Tofacitinib\nMA2025" ~ "Melon-Ardanaz 2025",
                  Cohort == "Infliximab\nT2011" ~ "Toedter 2011",
                )) %>%
  ggplot(aes(y = Source, x = Label)) +
  geom_point(aes(fill = NES, size = -log10(padj), alpha = Significance, shape = Direction)) +
  scale_shape_manual(values = direction_shapes) +
  scale_alpha_manual(values = significance_alpha) +
  scale_fill_gradient2(high = "coral3", low = "deepskyblue3", mid = "white") +
  facet_grid(Drug~Figure4, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

# C

cor_TOFvIFX <- cor.test(degs_R_W8vW0_TOFvIFX$log2FoldChange_TOF, degs_R_W8vW0_TOFvIFX$log2FoldChange_IFX)

scatterplot_degs_r_W8vW0_TOFvIFX_ggplotobj <- degs_R_W8vW0_TOFvIFX %>%
  ggplot(aes(x = log2FoldChange_TOF, y = log2FoldChange_IFX)) +
  geom_hline(yintercept = 0, color = "gray30") +
  geom_vline(xintercept = 0, color = "gray30") +
  geom_point_rast(size = 1.4, alpha = 0.8, col = "#808080") +
  geom_point_rast(data = . %>%
                    dplyr::filter(Significance_detailed == "NS"),
                  aes(color = Significance_detailed)) +
  geom_point_rast(data = . %>%
                    dplyr::filter(Significance_detailed != "NS"),
                  aes(color = Significance_detailed)) +
  geom_smooth(method = lm, level=0.99, col = "#000000") + 
  scale_color_manual(values = tofvifx_direction_colors,
                     name = "Significance") +
  annotate("text", 
           x = I(0.1), 
           y = I(0.95), 
           label = paste0("r = ", round(cor_TOFvIFX$estimate, 3))) + 
  labs(title = "Tofacitinib vs Infliximab",
       subtitle = "Responder: week 8 - week 0",
       x = "log2(w8-w0)\nTOF",
       y = "log2(w8-w0)\nIFX") +
  theme_bw() + 
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# D

## Own

rownames(vst_R_W8vW0_TOF) <- gsub("(ENSG[0-9]+)\\.[0-9]+$", "\\1", rownames(vst_R_W8vW0_TOF))

exprs_goi_R_W8vW0_TOF <- assay(vst_R_W8vW0_TOF)[rownames(vst_R_W8vW0_TOF) %in% goi_figure4$ENSG,]

sample_metadata_R_W8vW0_TOF <- colData(vst_R_W8vW0_TOF) %>%
  data.frame() %>%
  dplyr::mutate(DonorID = as.factor(DonorID)) %>%
  dplyr::select(Week, DonorID) %>%
  dplyr::arrange(Week, DonorID)

exprs_scaled_goi_R_W8vW0_TOF <- t(apply(as.matrix(exprs_goi_R_W8vW0_TOF), MARGIN = 1, FUN = scale))
colnames(exprs_scaled_goi_R_W8vW0_TOF) <- colnames(exprs_goi_R_W8vW0_TOF)
exprs_scaled_goi_R_W8vW0_TOF <- exprs_scaled_goi_R_W8vW0_TOF %>% 
  data.frame() %>%
  tibble::rownames_to_column(var = "ENSG") %>%
  dplyr::left_join(goi_figure4 %>%
                     dplyr::select(Gene, ENSG), by = "ENSG") %>%
  dplyr::select(-ENSG) %>%
  tibble::column_to_rownames(var = "Gene")

## Toedter 2011

sample_metadata_R_W8vW0_IFX <- colData(microarray_toedter2011_R_W8vW0_IFX_se) %>%
  data.frame() %>%
  dplyr::filter(SampleID != "GSM578741") %>%
  dplyr::select(Week, DonorID)

goi_w8mw0_R_W8vW0_IFX_toedter2011 <- assay(microarray_toedter2011_R_W8vW0_IFX_se, i = "logqnorm") %>%
  data.frame() %>%
  dplyr::mutate(gene_name = rowData(microarray_toedter2011_R_W8vW0_IFX_se)$`Gene Symbol`) %>%
  dplyr::filter(gene_name %in% goi_figure4$Gene) %>%
  tidyr::pivot_longer(-gene_name, names_to = "Sample_ID", values_to = "Expr") %>%
  dplyr::group_by(gene_name, Sample_ID) %>%
  dplyr::summarize(Expr = mean(Expr)) %>%
  tidyr::pivot_wider(id_cols = "gene_name", names_from = "Sample_ID", values_from = "Expr") %>%
  tibble::column_to_rownames(var = "gene_name")

goi_scaled_R_W8vW0_IFX_toedter2011 <- t(apply(as.matrix(goi_w8mw0_R_W8vW0_IFX_toedter2011), MARGIN = 1, FUN = scale))
colnames(goi_scaled_R_W8vW0_IFX_toedter2011) <- colnames(goi_w8mw0_R_W8vW0_IFX_toedter2011)

#### Combined

goi_scaled_R_W8vW0_IFXvTOF <- cbind(goi_scaled_R_W8vW0_IFX_toedter2011, exprs_scaled_goi_R_W8vW0_TOF[rownames(goi_scaled_R_W8vW0_IFX_toedter2011),])

sample_metadata_R_W8vW0_IFXvTOF <- sample_metadata_R_W8vW0_IFX %>% 
  dplyr::mutate(Drug = "IFX") %>%
  dplyr::rows_append(sample_metadata_R_W8vW0_TOF %>%
                       dplyr::mutate(Drug = "TOF")) %>%
  dplyr::arrange(Drug, Week)

feature_metadata_goi_R_W8vW0_TOFvIFX <- goi_figure4 %>%
  dplyr::filter(Gene %in% rownames(goi_scaled_R_W8vW0_IFXvTOF)) %>%
  tibble::column_to_rownames(var = "Gene")
feature_metadata_goi_R_W8vW0_TOFvIFX <- feature_metadata_goi_R_W8vW0_TOFvIFX[rownames(goi_scaled_R_W8vW0_IFXvTOF),]

cheatmap_sample_annotation_R_W8vW0_IFXvTOF <- rowAnnotation(
  Week = sample_metadata_R_W8vW0_IFXvTOF$Week,
  Drug = ifelse(sample_metadata_R_W8vW0_IFXvTOF$Drug == "IFX", "Infliximab", "Tofacitinib"),
  col = list(Week = timepoint_colors,
             Drug = treatment_colors),
  gp = gpar(col = "#D3D3D3")
)

cheatmap_goi_annotation_R_W8vW0_IFXvTOF <- HeatmapAnnotation(
  Geneset = feature_metadata_goi_R_W8vW0_TOFvIFX$Category, 
  col = list(Geneset = geneset_label_colors),
  gp = gpar(col = "#D3D3D3")
)

cheatmap_goi_exprs_scaled_R_W8vW0_IFXvTOF_ggplotobj <- ggplotify::as.ggplot(
  Heatmap(t(goi_scaled_R_W8vW0_IFXvTOF[,rownames(sample_metadata_R_W8vW0_IFXvTOF)]), 
        name = "Exprs", 
        row_split = sample_metadata_R_W8vW0_IFXvTOF$Drug,
        column_split = feature_metadata_goi_R_W8vW0_TOFvIFX$Category,
        top_annotation = cheatmap_goi_annotation_R_W8vW0_IFXvTOF,
        left_annotation = cheatmap_sample_annotation_R_W8vW0_IFXvTOF,
        cluster_rows = F,
        rect_gp = gpar(col = "#D3D3D3", lwd = 1),
        show_row_names = FALSE,
        row_title = NULL,
        column_title = NULL)
)

# Assembly

figA <- scatterplot_degs_r_W8vW0_ownvmelonardanaz2025_ggplotobj
figBD <- arrowplot_fgsea_RvNR_W8vW0_TOFvIFX_ownvmelonardanaz2025vtoedter2011_ggplotobj
figC <- scatterplot_degs_r_W8vW0_TOFvIFX_ggplotobj
figE <- cheatmap_goi_exprs_scaled_R_W8vW0_IFXvTOF_ggplotobj

# Save
ggsave(filename = figA_pdf, 
       plot = figA, 
       width = 6,
       height = 5.5)
ggsave(filename = figBD_pdf, 
       plot = figBD, 
       width = 6.5,
       height = 5.5)
ggsave(filename = figC_pdf, 
       plot = figC, 
       width = 8.75,
       height = 6)
ggsave(filename = figE_pdf, 
       plot = figE, 
       width = 10,
       height = 12.5)

sink(type = "message")
sink(type = "output")
close(log_file)