#!/usr/bin/env R
# Figure 5
# A) Scatterplot representing the longitudinal log2fold changes for colonic biopsies (X-axis) and peripheral blood (Y-axis) for responders.
# B) Arrowplot representing the pathways of interest over time for responders in PB and BINF.
# C) Heatmap visualization of genes representative of the NFKB canonical and non-canonical pathway for responders at week 0 and week 8.

# Paths
rnaseq_dds_rds <- snakemake@input[["rnaseq_dds_rds"]]
degs_BINF_R_W8vW0_csv <- snakemake@input[["degs_BINF_R_W8vW0_csv"]] #"output/260519/rna_Tissue_R_Inflamed_Week_8_vs_Baseline_deseq2_results.csv" 
degs_BINF_NR_W8vW0_csv <- snakemake@input[["degs_BINF_NR_W8vW0_csv"]] #"output/260519/rna_Tissue_NR_Inflamed_Week_8_vs_Baseline_deseq2_results.csv"
fgsea_BINF_R_W8vW0_csv <- snakemake@input[["fgsea_BINF_R_W8vW0_csv"]] #"output/260519/rna_Tissue_R_Inflamed_Week_8_vs_Baseline_fgsea_results.csv"
fgsea_BINF_NR_W8vW0_csv <- snakemake@input[["fgsea_BINF_NR_W8vW0_csv"]] #"output/260519/rna_Tissue_NR_Inflamed_Week_8_vs_Baseline_fgsea_results.csv"
degs_PB_R_W8vW0_csv <- snakemake@input[["degs_PB_R_W8vW0_csv"]] #"output/260519/rna_PBL_R_Inflamed_Week_8_vs_Baseline_deseq2_results.csv"
degs_PB_NR_W8vW0_csv <- snakemake@input[["degs_PB_NR_W8vW0_csv"]] #"output/260519/rna_PBL_NR_Inflamed_Week_8_vs_Baseline_deseq2_results.csv"
fgsea_PB_R_W8vW0_csv <- snakemake@input[["fgsea_PB_R_W8vW0_csv"]] #"output/260519/rna_PBL_R_Inflamed_Week_8_vs_Baseline_fgsea_results.csv"
fgsea_PB_NR_W8vW0_csv <- snakemake@input[["fgsea_PB_NR_W8vW0_csv"]] #"output/260519/rna_PBL_NR_Inflamed_Week_8_vs_Baseline_fgsea_results.csv"

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
library(DESeq2)
library(ggplot2)
library(scales)
library(ggrepel)
library(ggrastr)
library(patchwork)
library(ggplotify)
library(readxl)
library(readr)
library(scales)
library(stringr)
library(ComplexHeatmap)
library(ggh4x)
library(AnnotationDbi)
library(org.Hs.eg.db)

sessionInfo()

# Import
source(plotting_parameters_r)
goi <- readxl::read_excel(goi_xlsx)
pwoi <- readxl::read_excel(pwoi_xlsx)

dds <- readRDS(rnaseq_dds_rds)
print(degs_BINF_R_W8vW0_csv)
degs_BINF_R_W8vW0 <- read.csv(degs_BINF_R_W8vW0_csv)
print(degs_BINF_NR_W8vW0_csv)
degs_BINF_NR_W8vW0 <- read.csv(degs_BINF_NR_W8vW0_csv)
print(fgsea_BINF_R_W8vW0_csv)
fgsea_BINF_R_W8vW0 <- read.csv(fgsea_BINF_R_W8vW0_csv)
print(fgsea_BINF_NR_W8vW0_csv)
fgsea_BINF_NR_W8vW0 <- read.csv(fgsea_BINF_NR_W8vW0_csv)
print(degs_PB_R_W8vW0_csv)
degs_PB_R_W8vW0 <- read.csv(degs_PB_R_W8vW0_csv)
print(degs_PB_NR_W8vW0_csv)
degs_PB_NR_W8vW0 <- read.csv(degs_PB_NR_W8vW0_csv)
print(fgsea_PB_R_W8vW0_csv)
fgsea_PB_R_W8vW0 <- read.csv(fgsea_PB_R_W8vW0_csv)
print(fgsea_PB_NR_W8vW0_csv)
fgsea_PB_NR_W8vW0 <- read.csv(fgsea_PB_NR_W8vW0_csv)

# Preparation
goi_figure5 <- goi %>%
  dplyr::filter(!is.na(Figure5)) %>%
  dplyr::mutate(ENSG = mapIds(org.Hs.eg.db, 
                              keys = Gene,
                              column = "ENSEMBL", 
                              keytype = "SYMBOL"))
pwoi_figure5 <- pwoi %>%
  dplyr::filter(!is.na(Figure5))

## Combine degs PB R and NR
degs_PB_RvNR_W8vW0_combined <- degs_PB_R_W8vW0 %>%
  dplyr::inner_join(degs_PB_NR_W8vW0, by = c("ENSG", "ENSGv", "gene_symbol"), suffix = c("_R", "_NR")) %>%
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

## Combine degs BINF and PB R
degs_BINFvPB_R_W8vW0 <- degs_BINF_R_W8vW0 %>%
  dplyr::inner_join(degs_PB_R_W8vW0, 
                    by = c("ENSG", "ENSGv", "gene_symbol"), suffix = c("_BINF", "_PB")) %>%
  dplyr::mutate(Significance_detailed = factor(case_when(
    padj_BINF<0.05 & log2FoldChange_BINF<0 & pvalue_PB<0.05 & log2FoldChange_PB<0 ~ "Down",
    padj_BINF<0.05 & log2FoldChange_BINF>0 & pvalue_PB<0.05 & log2FoldChange_PB>0 ~ "Up",
    padj_BINF<0.05 & log2FoldChange_BINF<0 & pvalue_PB<0.05 & log2FoldChange_PB>0 ~ "B: down; PB: up",
    padj_BINF<0.05 & log2FoldChange_BINF>0 & pvalue_PB<0.05 & log2FoldChange_PB<0 ~ "B: up; PB: down",
    padj_BINF<0.05 & log2FoldChange_BINF<0 & pvalue_PB>0.05 ~ "B: down; PB: stable",
    padj_BINF<0.05 & log2FoldChange_BINF>0 & pvalue_PB>0.05 ~ "B: up; PB: stable",
    .default = "NS"
  ), levels = c("NS", "B: down; PB: up", "B: down; PB: stable", "B: up; PB: down", "B: up; PB: stable", "Down", "Up")))

degs_BINFvPB_R_W8vW0 %>% 
  dplyr::left_join(goi, by = c("gene_symbol" = "Gene")) %>%
  dplyr::filter(!is.na(Category)) %>%
  dplyr::select(gene_symbol, Category, log2FoldChange_BINF, log2FoldChange_PB, pvalue_BINF, pvalue_PB, Significance_detailed) %>%
  dplyr::arrange(Category, Significance_detailed)

degs_BINFvPB_R_W8vW0 %>% 
  dplyr::left_join(goi, by = c("gene_symbol" = "Gene")) %>%
  dplyr::filter(!is.na(Category)) %>%
  dplyr::select(gene_symbol, Category, log2FoldChange_BINF, log2FoldChange_PB, pvalue_BINF, pvalue_PB, Significance_detailed) %>%
  dplyr::arrange(Category, Significance_detailed) %>%
  dplyr::filter(Significance_detailed == "Down")

degs_BINFvPB_R_W8vW0 %>% 
  dplyr::left_join(goi, by = c("gene_symbol" = "Gene")) %>%
  dplyr::filter(!is.na(Category)) %>%
  dplyr::select(gene_symbol, Category, log2FoldChange_BINF, log2FoldChange_PB, pvalue_BINF, pvalue_PB, Significance_detailed) %>%
  dplyr::arrange(Category, Significance_detailed) %>%
  dplyr::filter(!Significance_detailed %in% c("B: down; PB: stable", "NS"))

# Combine degs BINF and PB R and NR for W8 and W0
fgsea_BINFvPB_RvNR_W8vW0 <- fgsea_BINF_R_W8vW0 %>%
  dplyr::select(pathway, pval, padj, NES) %>%
  dplyr::mutate(Tissue = "BINF",
                Response = "R") %>%
  dplyr::rows_append(fgsea_PB_R_W8vW0 %>%
                       dplyr::select(pathway, pval, padj, NES) %>%
                       dplyr::mutate(Tissue = "PB",
                                     Response = "R")) %>%
  dplyr::rows_append(fgsea_BINF_NR_W8vW0 %>%
                       dplyr::select(pathway, pval, padj, NES) %>%
                       dplyr::mutate(Tissue = "BINF",
                                     Response = "NR")) %>%
  dplyr::rows_append(fgsea_PB_NR_W8vW0 %>%
                       dplyr::select(pathway, pval, padj, NES) %>%
                       dplyr::mutate(Tissue = "PB",
                                     Response = "NR")) %>%
  dplyr::mutate(label = gsub("_", " ", pathway),
                label = gsub("(REACTOME|HALLMARK|KEGG)(.+)", "(\\1)\\2", label)) %>%
  dplyr::mutate(Direction = ifelse(NES<0, "Down", "Up"),
                Significance = ifelse(padj<0.05, "Significant", "NS"),
                Tissue = factor(Tissue, rev(c("BINF", "PB"))),
                Response = factor(Response, levels = c("NR", "R")))
  
# A

message("Starting on A")

cor_BINFvPB_W8vW0 <- cor.test(degs_BINFvPB_R_W8vW0$log2FoldChange_BINF, degs_BINFvPB_R_W8vW0$log2FoldChange_PB)

scatterplot_degs_binfvPB_r_W8vW0_ggplotobj <- degs_BINFvPB_R_W8vW0 %>%
  ggplot(aes(x = log2FoldChange_BINF, y = log2FoldChange_PB)) +
  geom_hline(yintercept = 0, color = "gray30") +
  geom_vline(xintercept = 0, color = "gray30") +
  geom_point_rast(data = . %>%
                    dplyr::filter(Significance_detailed == "NS"),
                  aes(color = Significance_detailed)) +
  geom_point_rast(data = . %>%
                    dplyr::filter(Significance_detailed != "NS"),
                  aes(color = Significance_detailed)) +
  scale_color_manual(values = bvpb_direction_colors,
                     name = "Significance") +
  scale_fill_manual(values = bvpb_direction_colors,
                    drop = FALSE,
                    guide = "none") +
  annotate("text", 
           x = I(0.1), 
           y = I(0.95), 
           label = paste0("r = ", round(cor_BINFvPB_W8vW0$estimate, 3))) + 
  labs(title = "Colonic biopsy vs peripheral blood",
       subtitle = "Responder: log2(week 8 - week 0)",
       x = "Colonic biopsy",
       y = "Peripheral blood") +
  theme_bw() + 
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# B

message("Starting on B")

arrowplot_fgsea_BINFvPB_R_W8vW0_ggplotobj <- fgsea_BINFvPB_RvNR_W8vW0 %>%
  dplyr::inner_join(pwoi_figure5, by = c("pathway"= "Geneset")) %>%
  dplyr::filter(Response == "R") %>%
  dplyr::mutate(Figure5 = factor(Figure5, levels = c("Down", "B: down"))) %>%
  ggplot(aes(y = Tissue, x = Label)) +
  geom_point(aes(fill = NES, size = -log10(padj), alpha = Significance, shape = Direction)) +
  scale_shape_manual(values = direction_shapes) +
  scale_alpha_manual(values = significance_alpha) +
  scale_fill_gradient2(high = "coral3", low = "deepskyblue3", mid = "white") +
  facet_grid(.~Figure5, scales = "free", space = "free") +
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

message("Starting on C")

PB_R_W8vW0_samples <- colData(dds) %>%
  data.frame() %>%
  dplyr::filter(Tissue == "PB",
                Response %in% c("R"),
                Week %in% c("W0", "W8"))

dds_PB_R_W8vW0 <- dds[,colnames(dds) %in% rownames(PB_R_W8vW0_samples)]
vst_PB_R_W8vW0 <- vst(dds_PB_R_W8vW0)
rownames(vst_PB_R_W8vW0) <- gsub("(ENSG[0-9]+)\\.[0-9]+$", "\\1", rownames(vst_PB_R_W8vW0))

goi_PB_R_W8vW0 <- assay(vst_PB_R_W8vW0)[goi_figure5$ENSG,]
goi_scaled_PB_R_W8vW0 <- t(apply(as.matrix(goi_PB_R_W8vW0), MARGIN = 1, FUN = scale))
colnames(goi_scaled_PB_R_W8vW0) <- colnames(goi_PB_R_W8vW0)

sample_metadata_PB_R_W8vW0 <- colData(vst_PB_R_W8vW0) %>%
  data.frame() %>%
  dplyr::arrange(Week, DonorID)

feature_metadata_goi <- data.frame(Geneset = goi_figure5$Category,
                                   Category = goi_figure5$Figure5,
                                   row.names = goi_figure5$ENSG)

cheatmap_goi_annotation_PB_R_W8vW0 <- HeatmapAnnotation(
  Geneset = feature_metadata_goi$Geneset,
  Category = feature_metadata_goi$Category,
  col = list(Geneset = geneset_label_colors,
             Category = category_label_colors),
  gp = gpar(col = "#D3D3D3")
)

cheatmap_sample_annotation_PB_R_W8vW0 <- rowAnnotation(
  Week = sample_metadata_PB_R_W8vW0$Week,
  col = list(Geneset = treatment_colors,
             Week = timepoint_colors),
  gp = gpar(col = "#D3D3D3")
)

cheatmap_PB_R_W8vW0_ggplotobj <- ggplotify::as.ggplot(
  Heatmap(t(goi_scaled_PB_R_W8vW0[,rownames(sample_metadata_PB_R_W8vW0)]), 
          row_split = sample_metadata_PB_R_W8vW0$Response,
          column_split = feature_metadata_goi$Geneset,
          top_annotation = cheatmap_goi_annotation_PB_R_W8vW0,
          left_annotation = cheatmap_sample_annotation_PB_R_W8vW0,
          rect_gp = gpar(col = "#D3D3D3", lwd = 1),
          column_labels = goi_figure5$Gene,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          name = "Expression",
          show_row_names = FALSE,
          row_title = NULL,
          column_title = NULL))

# Assembly
figA <- scatterplot_degs_binfvPB_r_W8vW0_ggplotobj
figB <- arrowplot_fgsea_BINFvPB_R_W8vW0_ggplotobj
figC <- cheatmap_PB_R_W8vW0_ggplotobj

# Save
ggsave(filename = figA_pdf, 
       plot = figA, 
       width = 6.5,
       height = 5.5)
ggsave(filename = figB_pdf, 
       plot = figB, 
       width = 7,
       height = 5.25)
ggsave(filename = figC_pdf, 
       plot = figC, 
       width = 9,
       height = 6)
