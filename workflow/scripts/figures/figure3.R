#!/usr/bin/env Rscript
# Figure 3
# A) Principal component analysis colored by response and week.
# B) Volcanoplots representing the longitudinal differentially expressed genes identified when comparing week 8 with week 0 in responders (top-right) and non-responders (bottom-left) as well as a scatterplot (bottom-right) representing the longitudinal log2fold changes for responders (X-axis) and non-responders (Y-axis).
# C) Scatterplot representing the longitudinal log2fold changes for responders (X-axis) and pretreatment differences between responders and non-responders pretreatment (Y-axis).
# D) Arrowplot representing the pathways of interest over time for responders and non-responders.
# E) Heatmap visualization of genes representative of pathways of interest over time for responders and non-responders.
# F) Principal component analysis of week 8 from responders, late responders, and non-responders.
# G) Heatmap visualization of genes representative of pathways of interest over time for late responders.
# H) Scatterplot of the log2(fold-change) calculated when comparing week 8 with week 0 for gene (X-axis) and protein (Y-axis) expression.
# I) Boxplot of the protein and gene expression of NFKB genes.

# Paths
rnaseq_dds_rds <- snakemake@input[["rnaseq_dds_rds"]] #"output/rnaseq/subsets/dds_general.rds"
degs_RvNR_W0_csv <- snakemake@input[["degs_RvNR_W0_csv"]] #"output/260519/rna_Tissue_Inflamed_Baseline_R_vs_NR_deseq2_results.csv" 
degs_R_W8vW0_csv <- snakemake@input[["degs_R_W8vW0_csv"]] #"output/260519/rna_Tissue_R_Inflamed_Week_8_vs_Baseline_deseq2_results.csv" 
degs_NR_W8vW0_csv <- snakemake@input[["degs_NR_W8vW0_csv"]] #"output/260519/rna_Tissue_NR_Inflamed_Week_8_vs_Baseline_deseq2_results.csv"
fgsea_RvNR_W0_csv <- snakemake@input[["fgsea_RvNR_W0_csv"]] #"output/260519/rna_Tissue_Inflamed_Baseline_R_vs_NR_fgsea_results.csv"
fgsea_R_W8vW0_csv <- snakemake@input[["fgsea_R_W8vW0_csv"]] #"output/260519/rna_Tissue_R_Inflamed_Week_8_vs_Baseline_fgsea_results.csv"
fgsea_NR_W8vW0_csv <- snakemake@input[["fgsea_NR_W8vW0_csv"]] #"output/260519/rna_Tissue_NR_Inflamed_Week_8_vs_Baseline_fgsea_results.csv"

olink_se_rds <- snakemake@input[["olink_se_rds"]] #"output/olink/exprs.csv"
deps_R_W8vW0_csv <- snakemake@input[["deps_R_W8vW0_csv"]] #"output/260519/olink_Tissue_R_Inflamed_Week_8_vs_Baseline_limma_results.csv"
deps_NR_W8vW0_csv <- snakemake@input[["deps_NR_W8vW0_csv"]] #"output/260519/olink_Tissue_NR_Inflamed_Week_8_vs_Baseline_limma_results.csv"

goi_xlsx <- snakemake@params[["goi_xlsx"]] #"config/features_of_interest/260716_genes.xlsx"
pwoi_xlsx <- snakemake@params[["pwoi_xlsx"]] #"config/features_of_interest/260716_genesets.xlsx"
plotting_parameters_r <- snakemake@params[["plotting_parameters"]] #"workflow/scripts/plotting_parameters.R"

figA_pdf <- snakemake@output[["figA_pdf"]]
figB_pdf <- snakemake@output[["figB_pdf"]]
figC_pdf <- snakemake@output[["figC_pdf"]]
figD_pdf <- snakemake@output[["figD_pdf"]]
figEG_pdf <- snakemake@output[["figEG_pdf"]]
figF_pdf <- snakemake@output[["figF_pdf"]]
figH_pdf <- snakemake@output[["figH_pdf"]]
figI_pdf <- snakemake@output[["figI_pdf"]]

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
print(degs_RvNR_W0_csv)
degs_RvNR_W0 <- read.csv(degs_RvNR_W0_csv)
print(degs_R_W8vW0_csv)
degs_R_W8vW0 <- read.csv(degs_R_W8vW0_csv)
print(degs_NR_W8vW0_csv)
degs_NR_W8vW0 <- read.csv(degs_NR_W8vW0_csv)
print(fgsea_RvNR_W0_csv)
fgsea_RvNR_W0 <- read.csv(fgsea_RvNR_W0_csv)
print(fgsea_R_W8vW0_csv)
fgsea_R_W8vW0 <- read.csv(fgsea_R_W8vW0_csv)
print(fgsea_NR_W8vW0_csv)
fgsea_NR_W8vW0 <- read.csv(fgsea_NR_W8vW0_csv)

olink_se <- readRDS(olink_se_rds)
print(deps_R_W8vW0_csv)
deps_R_W8vW0 <- read.csv(deps_R_W8vW0_csv)
print(deps_NR_W8vW0_csv)
deps_NR_W8vW0 <- read.csv(deps_NR_W8vW0_csv)

# Preparation

goi_figure3 <- goi %>%
  dplyr::filter(!is.na(Figure3)) %>%
  dplyr::mutate(ENSG = mapIds(org.Hs.eg.db, 
                              keys = Gene,
                              column = "ENSEMBL", 
                              keytype = "SYMBOL"))
pwoi_figure3 <- pwoi %>%
  dplyr::filter(!is.na(Figure3))

## Significant results RvNR_W0

sig_degs_RvNR_W0 <- degs_RvNR_W0 %>%
  dplyr::filter(padj<0.05) %>%
  dplyr::arrange(pvalue)

sig_fgsea_RvNR_W0 <- fgsea_RvNR_W0 %>%
  dplyr::filter(padj<0.05) %>%
  dplyr::arrange(pval)

## Combine degs R and NR
degs_RvNR_W8vW0_combined <- degs_R_W8vW0 %>%
  dplyr::inner_join(degs_NR_W8vW0, by = c("ENSG", "ENSGv", "gene_symbol"), suffix = c("_R", "_NR")) %>%
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

## Combine degs and deps R
depvdeg_R_W8vW0_combined <- degs_R_W8vW0 %>%
  dplyr::inner_join(deps_R_W8vW0 %>%
                      dplyr::group_by(gene_symbol) %>%
                      dplyr::summarize(log2FoldChange = mean(log2FoldChange),
                                       AveExpr = mean(AveExpr),
                                       stat = mean(stat),
                                       pvalue = mean(pvalue)) %>%
                      dplyr::mutate(padj = p.adjust(pvalue, method = "BH")),
                    by = c("gene_symbol"), 
                    suffix = c("_GEX", "_PEX")) %>%
  dplyr::mutate(Significance_detailed = case_when(
    padj_GEX<0.05 & log2FoldChange_GEX<0 & pvalue_PEX<0.05 & log2FoldChange_PEX<0 ~ "GEX & PEX down",
    padj_GEX<0.05 & log2FoldChange_GEX>0 & pvalue_PEX<0.05 & log2FoldChange_PEX>0 ~ "GEX & PEX up",
    .default = "NS"
  ),
  Significance = case_when(
    Significance_detailed == "GEX & PEX down" ~ "Down",
    Significance_detailed == "GEX & PEX up" ~ "Up",
    .default = "NS"
  ),
  Significance = factor(Significance, levels = names(direction_colors)))

## Combine degs RvNR and R:W8vW0 
degs_RvNRW0_RW8vW0_combined <- degs_RvNR_W0 %>%
  dplyr::inner_join(degs_R_W8vW0,
                    by = c("ENSG", "ENSGv", "gene_symbol"), suffix = c("_W0RvNR", "_RW8vW0")) %>%
  dplyr::mutate(Significance = case_when(
    padj_W0RvNR<0.05 & log2FoldChange_W0RvNR<0 & pvalue_RW8vW0<0.05 & log2FoldChange_RW8vW0<0 ~ "W0: low; W0>W8: down",
    padj_W0RvNR<0.05 & log2FoldChange_W0RvNR>0 & pvalue_RW8vW0<0.05 & log2FoldChange_RW8vW0<0 ~ "W0: high; W0>W8: down",
    padj_W0RvNR<0.05 & log2FoldChange_W0RvNR<0 & pvalue_RW8vW0<0.05 & log2FoldChange_RW8vW0>0 ~ "W0: low; W0>W8: up",
    padj_W0RvNR<0.05 & log2FoldChange_W0RvNR>0 & pvalue_RW8vW0<0.05 & log2FoldChange_RW8vW0>0 ~ "W0: high; W0>W8: up",
    padj_W0RvNR<0.05 & log2FoldChange_W0RvNR<0 & pvalue_RW8vW0>0.05 ~ "W0: low; W0>W8: stable",
    padj_W0RvNR<0.05 & log2FoldChange_W0RvNR>0 & pvalue_RW8vW0>0.05 ~ "W0: high; W0>W8: stable",
    .default = "NS"
  ))

## Combine degs R and NR
fgsea_RvNR_W8vW0 <- fgsea_R_W8vW0 %>%
  dplyr::mutate(Response = "R") %>%
  dplyr::rows_append(fgsea_NR_W8vW0 %>%
                       dplyr::mutate(Response = "NR")) %>%
  dplyr::mutate(Direction = ifelse(NES<0, "Down", "Up"),
                Significance = ifelse(padj<0.05, "Significant", "NS"),
                Response = factor(Response, levels = c("R", "NR")))

# A

message("Starting on A")

## Filter samples of interest
samples_RvNR_W8vW0 <- colData(dds) %>%
  data.frame() %>%
  dplyr::filter(Tissue == "BINF",
                Response %in% c("R", "NR"),
                Week %in% c("W0", "W8"))

dds_RvNR_W8vW0 <- dds[,which(colnames(dds) %in% rownames(samples_RvNR_W8vW0))]
dds_RvNR_W8vW0 <- dds_RvNR_W8vW0[which(rowSums(counts(dds_RvNR_W8vW0) >= 10) >= 5),]
vst_RvNR_W8vW0 <- vst(dds_RvNR_W8vW0)
rownames(vst_RvNR_W8vW0) <- gsub("(ENSG[0-9]+)\\.[0-9]+$", "\\1", rownames(vst_RvNR_W8vW0))

svd_vst_RvNR_W8vW0 <- svd(t(assay(vst_RvNR_W8vW0)))
svd_vst_RvNR_W8vW0_df <- data.frame(PC1 = svd_vst_RvNR_W8vW0$u[,1],
                                    PC2 = svd_vst_RvNR_W8vW0$u[,2],
                                    Response = colData(vst_RvNR_W8vW0)$Response,
                                    Week = colData(vst_RvNR_W8vW0)$Week,
                                    DonorID = colData(vst_RvNR_W8vW0)$DonorID,
                                    Sex = colData(vst_RvNR_W8vW0)$Sex,
                                    row.names = colnames(vst_RvNR_W8vW0))

pca_RvNR_W8vW0_split_ggplotobj <- ggplot(svd_vst_RvNR_W8vW0_df, aes(x = PC1, y = PC2)) +
  stat_ellipse(alpha = 0.5, linetype = "dashed", aes(col = Response)) +
  geom_point(shape = 21, size = 3, aes(fill = Response)) +
  facet_wrap(~Week) +
  scale_color_manual(values = response_colors) +
  scale_fill_manual(values = response_colors) +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

pca_RvNR_W8vW0_ggplotobj <- ggplot(svd_vst_RvNR_W8vW0_df, aes(x = PC1, y = PC2)) +
  stat_ellipse(alpha = 0.5, linetype = "dashed", aes(col = Response)) +
  geom_point(shape = 21, size = 3, aes(fill = Response, alpha = Week)) +
  scale_color_manual(values = response_colors) +
  scale_fill_manual(values = response_colors) +
  scale_alpha_manual(values = c("W0" = 0.25, "W8" = 1)) +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# B

message("Starting on B")

## Volcanoplot R:W8vW0

l2fc_R_W8vW0_range <- range(degs_RvNR_W8vW0_combined$log2FoldChange_R)

volcanoplot_degs_R_W8vW0_ggplotobj <- degs_RvNR_W8vW0_combined %>%
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
       y = "-log10(padj)") +
  guides(fill=guide_legend(title="R")) +
  scale_y_continuous(position = "right") +
  coord_cartesian(xlim = l2fc_R_W8vW0_range) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

## Volcanoplot NR:W8vW0

l2fc_NR_W8vW0_range <- range(degs_RvNR_W8vW0_combined$log2FoldChange_NR)

volcanoplot_degs_NR_W8vW0_ggplotobj <- degs_RvNR_W8vW0_combined %>%
  dplyr::mutate(mlog10padj = -log10(padj_NR),
                Significance = case_when(
                  log2FoldChange_NR > 0 & padj_NR < 0.05 ~ "Up",
                  log2FoldChange_NR < 0 & padj_NR < 0.05 ~ "Down",
                  .default = "NS")) %>%
  ggplot(aes(y = log2FoldChange_NR, x = mlog10padj, color = Significance)) +
  scale_x_reverse() +
  geom_vline(xintercept = 0, color = "gray40") +
  geom_hline(yintercept = 0, color = "gray40") +
  geom_point_rast(size = 1.4, alpha = 0.8) +
  scale_color_manual(values = direction_colors,
                     drop = FALSE) +
  labs(y = "log2(fold change)", 
       x = "-log10(padj)") +
  coord_cartesian(ylim = l2fc_NR_W8vW0_range) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

## Scatterplot RvNR:W8vW0

cor_RvNR <- cor.test(degs_RvNR_W8vW0_combined$log2FoldChange_R, degs_RvNR_W8vW0_combined$log2FoldChange_NR)

scatterplot_degs_RvNR_W8vW0_ggplotobj <- degs_RvNR_W8vW0_combined %>%
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
  scale_y_continuous(position = "right") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

degs_RvNR_W8vW0_ggplotobj <- guide_area() + 
  volcanoplot_degs_R_W8vW0_ggplotobj +
  volcanoplot_degs_NR_W8vW0_ggplotobj + 
  scatterplot_degs_RvNR_W8vW0_ggplotobj  +
  plot_layout(guides = "collect",
              widths = c(1, 1), 
              heights = c(1, 1))

# C

message("Starting on C")

fgsea_RvNR_W8vW0_ggplotobj <- fgsea_RvNR_W8vW0 %>%
  dplyr::inner_join(pwoi_figure3 %>% 
                      dplyr::rename(Group = Figure3), 
                    by = c("pathway" = "Geneset")) %>%
  dplyr::mutate(Group = factor(Group, levels = c("R: down", "R: up", "R: down; NR: up", "Down"))) %>%
  ggplot(aes(x = forcats::fct_reorder(Label, -log10(pval)), y = Response)) +
  geom_point(aes(fill = NES, size = -log10(padj), alpha = Significance, shape = Direction)) +
  scale_shape_manual(values = direction_shapes) +
  scale_alpha_manual(values = significance_alpha) +
  scale_fill_gradient2(high = "coral3", low = "deepskyblue3", mid = "white") +
  scale_y_discrete(labels = label_wrap(50)) +
  facet_grid(.~Group, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

# Dv1

message("Starting on D")

goi_RvNR_W8vW0 <- assay(vst_RvNR_W8vW0)[goi_figure3$ENSG,]
goi_scaled_RvNR_W8vW0 <- t(apply(as.matrix(goi_RvNR_W8vW0), MARGIN = 1, FUN = scale))
colnames(goi_scaled_RvNR_W8vW0) <- colnames(goi_RvNR_W8vW0)

sample_metadata_RvNR_W8vW0 <- colData(vst_RvNR_W8vW0) %>%
  data.frame() %>%
  dplyr::mutate(Response_timepoint = paste0(Response, Week))

feature_metadata_goi <- data.frame(Category = goi_figure3$Category,
                                                 row.names = goi_figure3$ENSG)

cheatmap_goi_annotation_RvNR_W8vW0 <- HeatmapAnnotation(
  Geneset = feature_metadata_goi$Category, 
  col = list(Geneset = geneset_label_colors),
  gp = gpar(col = "#D3D3D3")
)

cheatmap_sample_annotation_RvNR_W8vW0 <- rowAnnotation(
  Week = sample_metadata_RvNR_W8vW0$Week,
  Response = sample_metadata_RvNR_W8vW0$Response,
  col = list(Geneset = treatment_colors,
             Response = response_colors,
             Week = timepoint_colors),
  gp = gpar(col = "#D3D3D3")
)

cheatmap_RvNR_W8vW0 <- Heatmap(t(goi_scaled_RvNR_W8vW0[,rownames(sample_metadata_RvNR_W8vW0)]), 
        row_split = sample_metadata_RvNR_W8vW0$Response,
        column_split = feature_metadata_goi$Category,
        top_annotation = cheatmap_goi_annotation_RvNR_W8vW0,
        left_annotation = cheatmap_sample_annotation_RvNR_W8vW0,
        rect_gp = gpar(col = "#D3D3D3", lwd = 1),
        column_labels = goi_figure3$Gene,
        name = "Expression",
        show_row_names = FALSE,
        row_title = NULL,
        column_title = NULL)

# E

message("Starting on E")

## Filter samples of interest
samples_RvLRvNR_W8 <- colData(dds) %>%
  data.frame() %>%
  dplyr::filter(Tissue == "BINF",
                Response %in% c("R", "LR", "NR"),
                Week %in% c("W8"))

# Filter genes with at least 10 counts in less than 5 samples
dds_RvLRvNR_W8 <- dds[,colnames(dds) %in% rownames(samples_RvLRvNR_W8)]
dds_RvLRvNR_W8 <- dds_RvLRvNR_W8[rowSums(counts(dds_RvLRvNR_W8) >= 10) >= 5,]
vst_RvLRvNR_W8 <- vst(dds_RvLRvNR_W8)

svd_vst_RvLRvNR_W8 <- svd(t(assay(vst_RvLRvNR_W8)))
svd_vst_RvLRvNR_W8_df <- data.frame(PC1 = svd_vst_RvLRvNR_W8$u[,1],
                                    PC2 = svd_vst_RvLRvNR_W8$u[,2],
                                    Response = factor(colData(vst_RvLRvNR_W8)$Response, levels = c("R", "LR", "NR")), 
                                    Week = colData(vst_RvLRvNR_W8)$Week,
                                    DonorID = colData(vst_RvLRvNR_W8)$DonorID,
                                    Sex = colData(vst_RvLRvNR_W8)$Sex,
                                    row.names = colnames(vst_RvLRvNR_W8))

pca_RvLRvNR_W8_ggplotobj <- svd_vst_RvLRvNR_W8_df %>% 
  dplyr::filter(Week == "W8") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  stat_ellipse(linetype = "dashed", aes(col = Response, alpha = Response)) +
  geom_point(shape = 21, size = 3, aes(fill = Response, alpha = Response)) +
  # facet_wrap(~Week) +
  scale_color_manual(values = response_colors) +
  scale_fill_manual(values = response_colors) +
  scale_alpha_manual(values = c("R" = 0.4, "LR" = 1, "NR" = 0.4)) +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# Dv3

cor_RvNRW0_RW8vW0 <- cor.test(degs_RvNRW0_RW8vW0_combined$log2FoldChange_W0RvNR, degs_RvNRW0_RW8vW0_combined$log2FoldChange_RW8vW0)

scatterplot_degs_RvNRW0_RW8vW0_ggplotobj <- degs_RvNRW0_RW8vW0_combined %>%
  ggplot(aes(x = log2FoldChange_RW8vW0, y = log2FoldChange_W0RvNR)) +
  geom_hline(yintercept = 0, color = "gray40") +
  geom_vline(xintercept = 0, color = "gray40") +
  geom_point_rast(aes(color = Significance), col = "#808080") +
  geom_point_rast(data = . %>%
                    dplyr::filter(Significance == "NS"),
                  col = "#d3d3d3") +
  geom_point_rast(data = . %>%
                    dplyr::filter(Significance != "NS"),
                  col = "#000000") +
  geom_label_repel(data = . %>%
                     dplyr::filter(gene_symbol %in% (goi %>% 
                                                       dplyr::filter(!is.na(Figure2)) %>%
                                                       dplyr::pull(Gene))),
                   aes(label = gene_symbol),
                   size = 2.8,
                   max.overlaps = Inf,
                   force = 3,
                   force_pull = 0.5,
                   max.iter = 10000,
                   segment.colour = "gray50",
                   segment.size = 0.6,
                   segment.alpha = 0.9,
                   point.padding = unit(0.6, "lines"),
                   box.padding = unit(0.4, "lines"),
                   min.segment.length = 0,
                   fontface = "bold", # Bold for contrast
                   label.padding = unit(0.15, "lines"),
                   show.legend = FALSE) +
  labs(x = "R: W8-W0",
       y = "W0: R-NR") +
  annotate("text", 
           x = I(0.1), 
           y = I(0.95), 
           label = paste0("r = ", round(cor_RvNRW0_RW8vW0$estimate, 3))) + 
  theme_bw() + 
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# D/F 

message("Starting on D and F")

## Filter samples of interest
samples_RvLRvNR_W8vW0 <- colData(dds) %>%
  data.frame() %>%
  dplyr::filter(Tissue == "BINF",
                Response %in% c("R", "LR", "NR"),
                Week %in% c("W0", "W8"))

# Filter genes with at least 10 counts in less than 5 samples
dds_RvLRvNR_W8vW0 <- dds[,which(colnames(dds) %in% rownames(samples_RvLRvNR_W8vW0))]
dds_RvLRvNR_W8vW0 <- dds_RvLRvNR_W8vW0[which(rowSums(counts(dds_RvLRvNR_W8vW0) >= 10) >= 5),]
vst_RvLRvNR_W8vW0 <- vst(dds_RvLRvNR_W8vW0)
rownames(vst_RvLRvNR_W8vW0) <- gsub("(ENSG[0-9]+)\\.[0-9]+$", "\\1", rownames(vst_RvLRvNR_W8vW0))

goi_RvLRvNR_W8vW0 <- assay(vst_RvLRvNR_W8vW0)[goi_figure3$ENSG,]
goi_scaled_RvLRvNR_W8vW0 <- t(apply(as.matrix(goi_RvLRvNR_W8vW0), MARGIN = 1, FUN = scale))
colnames(goi_scaled_RvLRvNR_W8vW0) <- colnames(goi_RvLRvNR_W8vW0)

sample_metadata_RvLRvNR_W8vW0 <- colData(vst_RvLRvNR_W8vW0) %>%
  data.frame() %>%
  dplyr::mutate(Response_timepoint = paste0(Response, Week))

feature_metadata_goi <- data.frame(Category = goi_figure3$Category,
                                   row.names = goi_figure3$ENSG)

cheatmap_goi_annotation_RvLRvNR_W8vW0 <- HeatmapAnnotation(
  Geneset = feature_metadata_goi$Category, 
  col = list(Geneset = geneset_label_colors),
  gp = gpar(col = "#D3D3D3")
)

cheatmap_sample_annotation_RvLRvNR_W8vW0 <- rowAnnotation(
  Week = sample_metadata_RvLRvNR_W8vW0$Week,
  Response = sample_metadata_RvLRvNR_W8vW0$Response,
  col = list(Geneset = treatment_colors,
             Response = response_colors,
             Week = timepoint_colors),
  gp = gpar(col = "#D3D3D3")
)

cheatmap_goi_scaled_RvLRvNR_W8vW0_ggplotobj <- ggplotify::as.ggplot(
  Heatmap(t(goi_scaled_RvLRvNR_W8vW0[,rownames(sample_metadata_RvLRvNR_W8vW0)]), 
          row_split = factor(sample_metadata_RvLRvNR_W8vW0$Response, levels = c("LR", "R", "NR")), #OK
          # row_split = factor(sample_metadata_RvLRvNR_W8vW0$Response, levels = c("R", "LR", "NR")), #NO
          # row_split = factor(sample_metadata_RvLRvNR_W8vW0$Response, levels = c("NR", "R", "LR")), #NO
          # row_split = factor(sample_metadata_RvLRvNR_W8vW0$Response, levels = c("R", "NR", "LR")), #NO
          column_split = feature_metadata_goi$Category,
          top_annotation = cheatmap_goi_annotation_RvLRvNR_W8vW0,
          left_annotation = cheatmap_sample_annotation_RvLRvNR_W8vW0,
          rect_gp = gpar(col = "#D3D3D3", lwd = 1),
          column_labels = goi_figure3$Gene,
          name = "Expression",
          show_row_names = FALSE,
          row_title = NULL,
          column_title = NULL))

# G

message("Starting on G")

cor_depvdeg <- cor.test(depvdeg_R_W8vW0_combined$log2FoldChange_GEX, depvdeg_R_W8vW0_combined$log2FoldChange_PEX)

scatterplot_depvdeg_R_W8vW0_ggplotobj <- depvdeg_R_W8vW0_combined %>%
  ggplot(aes(x = log2FoldChange_GEX, y = log2FoldChange_PEX)) +
  geom_hline(yintercept = 0, color = "gray40") +
  geom_vline(xintercept = 0, color = "gray40") +
  geom_point_rast(data = . %>%
                    dplyr::filter(Significance == "NS"),
                  aes(color = Significance)) +
  geom_point_rast(data = . %>%
                    dplyr::filter(Significance != "NS"),
                  aes(color = Significance)) +
  geom_label_repel(data = . %>%
                     dplyr::filter(gene_symbol %in% goi_figure3$Gene),
                   aes(label = gene_symbol, col = Significance),  
                   nudge_y = 0.1,
                   size = 2.8,
                   max.overlaps = Inf,
                   force = 3, 
                   force_pull = 0.5,
                   max.iter = 10000,
                   segment.colour = "gray50",
                   segment.size = 0.6,
                   segment.alpha = 0.9,
                   point.padding = unit(0.6, "lines"),
                   box.padding = unit(0.4, "lines"),
                   min.segment.length = 0,
                   fontface = "bold",
                   label.padding = unit(0.15, "lines"),
                   show.legend = FALSE) +
  scale_color_manual(values = direction_colors, 
                     name = "Significance") +
  scale_fill_manual(values = direction_colors,
                    drop = FALSE,
                    guide = "none") +
  labs(x = "log2(w8-w0)\nRNA-seq",
       y = "log2(w8-w0)\nOlink") +
  annotate("text", 
           x = I(0.1), 
           y = I(0.95), 
           label = paste0("r = ", round(cor_depvdeg$estimate, 3))) + 
  theme_bw() + 
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# H

message("Starting on H")

selected_features_olink_rnaseq <- c("IL6", "ICAM1", "CXCL8", "IL1B", "TNFRSF8")

samples_olink_RvNR_W8vW0 <- colData(olink_se) %>%
  data.frame() %>%
  dplyr::filter(Tissue == "BINF",
                Response %in% c("R", "NR"),
                Week %in% c("W0", "W8"))

olink_RvNR_W8vW0_se <- olink_se[,rownames(samples_olink_RvNR_W8vW0)]

olink_selected_features <- rowData(olink_RvNR_W8vW0_se) %>%
  data.frame() %>%
  dplyr::filter(Assay %in% selected_features_olink_rnaseq)

exprs_olink_RvNR_W8vW0 <- assay(olink_RvNR_W8vW0_se)[rownames(olink_selected_features),] %>%
  tibble::rownames_to_column(var = "Feature_ID") %>%
  tidyr::pivot_longer(-Feature_ID, values_to = "Exprs", names_to = "Sample_ID") %>%
  dplyr::left_join(colData(olink_RvNR_W8vW0_se) %>%
                     data.frame() %>%
                     tibble::rownames_to_column(var = "Sample_ID"),
                   by = "Sample_ID") %>%
  dplyr::left_join(rowData(olink_RvNR_W8vW0_se) %>%
                     data.frame() %>%
                     tibble::rownames_to_column(var = "Feature_ID")) %>%
  dplyr::select(Feature_ID, Sample_ID, DonorID, Response, Week, Assay, Exprs) %>%
  dplyr::mutate(Modality = "Olink") %>%
  dplyr::group_by(Sample_ID, DonorID, Response, Week, Assay, Modality) %>%
  dplyr::summarize(Exprs = mean(Exprs)) %>%
  dplyr::rename(Feature = Assay)

rnaseq_selected_features <- goi_figure3 %>% 
  dplyr::filter(Gene %in% selected_features_olink_rnaseq)

exprs_rnaseq_RvNR_W8vW0 <- assay(vst_RvNR_W8vW0)[rnaseq_selected_features$ENSG,] %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "Feature_ID") %>%
  tidyr::pivot_longer(-Feature_ID, values_to = "Exprs", names_to = "Sample_ID") %>%
  dplyr::left_join(colData(vst_RvNR_W8vW0) %>%
                     data.frame() %>%
                     tibble::rownames_to_column(var = "Sample_ID"),
                   by = "Sample_ID") %>%
  dplyr::left_join(goi_figure3, 
                   by = c("Feature_ID" = "ENSG")) %>%
  dplyr::mutate(Modality = "RNA-seq") %>%
  dplyr::select(Sample_ID, DonorID, Response, Week, Gene, Modality, Exprs) %>%
  dplyr::rename(Feature = Gene)

exprs_olinknrnaseq_RvNR_W8vW0 <- exprs_rnaseq_RvNR_W8vW0 %>%
  rows_append(exprs_olink_RvNR_W8vW0) %>%
  dplyr::mutate(Modality = factor(Modality, levels = c("RNA-seq", "Olink")))

boxplot_goi_olinknrnaseq_RvNR_W8vW0_ggplotobj <- exprs_olinknrnaseq_RvNR_W8vW0 %>%
  ggplot(aes(x = Week, y = Exprs)) +
  ggh4x::facet_grid2(Modality~Feature, independent = "y", scales = "free_y") +
  geom_boxplot(aes(fill = Response), alpha = 0.25, outlier.shape = NA) +
  geom_point(aes(group = DonorID, shape = Response, fill = Response), show.legend = F, position = position_dodge(0.1)) +
  geom_line(aes(group = DonorID, alpha = 0.75, col = Response), linetype = "dotted", show.legend = F, position = position_dodge(0.1)) +
  labs(x = "Timepoint",
       y = "Expression") +
  scale_color_manual(values = response_colors) + 
  scale_fill_manual(values = response_colors) + 
  scale_shape_manual(values = response_shapes) +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

# Assembly

figA <- pca_RvNR_W8vW0_ggplotobj
figB <- degs_RvNR_W8vW0_ggplotobj
figC <- scatterplot_degs_RvNRW0_RW8vW0_ggplotobj
figD <- fgsea_RvNR_W8vW0_ggplotobj
figEG <- cheatmap_goi_scaled_RvLRvNR_W8vW0_ggplotobj
figF <- pca_RvLRvNR_W8_ggplotobj
figH <- scatterplot_depvdeg_R_W8vW0_ggplotobj
figI <- boxplot_goi_olinknrnaseq_RvNR_W8vW0_ggplotobj

# Save

ggsave(filename = figA_pdf, 
       plot = figA, 
       width = 5,
       height = 5.5)
ggsave(filename = figB_pdf, 
       plot = figB, 
       width = 6,
       height = 6.5)
ggsave(filename = figC_pdf,
       plot = figC,
       width = 6,
       height = 6.5)
ggsave(filename = figD_pdf, 
       plot = figD, 
       width = 6.5,
       height = 5.25)
ggsave(filename = figEG_pdf, 
       plot = figEG, 
       width = 11,
       height = 12.5)
ggsave(filename = figF_pdf, 
       plot = figF, 
       width = 5,
       height = 5.5)
ggsave(filename = figH_pdf, 
       plot = figH, 
       width = 6,
       height = 5)
ggsave(filename = figI_pdf, 
       plot = figI, 
       width = 8,
       height = 4)

sink(type = "message")
sink(type = "output")
close(log_file)