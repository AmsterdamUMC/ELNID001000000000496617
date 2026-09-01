#!/usr/bin/env Rscript
# Figure 2
# A) Principal component analysis of future responders and non-responders pretreatment.
# B) Volcanoplots representing the longitudinal differentially expressed genes identified when comparing responders with non-responders pretreatment.
# C) Heatmap visualization of genes STING1, IL36G, IL36RN, PRKCQ, and NFATC2

# Paths
sample_metadata_xlsx <- snakemake@input[["sample_metadata_xlsx"]] #"config/metadata/samples/Sample_metadata_all.xlsx"
dds_rds <- snakemake@input[["dds_rds"]] #"output/rnaseq/subsets/dds_general.rds"
degs_csv <- snakemake@input[["degs_csv"]] #"output/260519/rna_Tissue_Inflamed_Baseline_R_vs_NR_deseq2_results.csv" #args[1]
fgsea_csv <- snakemake@input[["fgsea_csv"]] #"output/260519/rna_Tissue_Inflamed_Baseline_R_vs_NR_fgsea_results.csv" #args[1]
goi_xlsx <- snakemake@params[["goi_xlsx"]] #"config/features_of_interest/260716_genes.xlsx"
pwoi_xlsx <- snakemake@params[["pwoi_xlsx"]] #"config/features_of_interest/260716_genesets.xlsx"
plotting_parameters_r <- snakemake@params[["plotting_parameters"]] #"workflow/scripts/plotting_parameters.R"
figA_pdf <- snakemake@output[["figA_pdf"]]
figB_pdf <- snakemake@output[["figB_pdf"]]
figC_pdf <- snakemake@output[["figC_pdf"]]
figD_pdf <- snakemake@output[["figD_pdf"]]
log_file <- snakemake@log[[1]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

print(sample_metadata_xlsx)
print(dds_rds)
print(degs_csv)
print(fgsea_csv)
print(goi_xlsx)
print(pwoi_xlsx)
print(plotting_parameters_r)

# Libraries
library(dplyr)
library(DESeq2)
library(ggplot2)
library(scales)
library(ggrepel)
library(ggrastr)
library(patchwork)
library(readxl)
library(readr)
library(stringr)
library(ComplexHeatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(MASS)

sessionInfo()

# Import
source(plotting_parameters_r)

sample_metadata <- readxl::read_excel(sample_metadata_xlsx)
goi <- readxl::read_excel(goi_xlsx)
pwoi <- readxl::read_excel(pwoi_xlsx)
dds <- readRDS(dds_rds)
degs_BINF_RvNR_W0_TOF <- read.csv(degs_csv)
fgsea_BINF_RvNR_W0_TOF <- read.csv(fgsea_csv)

# Prepare features of interest
goi_figure2 <- goi %>%
  dplyr::filter(!is.na(Figure2)) %>%
  dplyr::mutate(ENSG = mapIds(org.Hs.eg.db, 
                              keys = Gene,
                              column = "ENSEMBL", 
                              keytype = "SYMBOL"))
pwoi_figure2 <- pwoi %>%
  dplyr::filter(!is.na(Figure2))

sig_degs_BINF_RvNR_W0_TOF <- degs_BINF_RvNR_W0_TOF %>%
  dplyr::filter(padj<0.05) %>%
  dplyr::arrange(pvalue)

sig_fgsea_BINF_RvNR_W0_TOF <- fgsea_BINF_RvNR_W0_TOF %>%
  dplyr::filter(padj<0.05) %>%
  dplyr::arrange(pval)

# A

## Filter samples of interest
BINF_RvNR_W0_TOF_samples <- colData(dds) %>%
  data.frame() %>%
  dplyr::filter(Tissue == "BINF",
                Response %in% c("R", "NR"),
                Week %in% c("W0"))

# Filter genes with at least 10 counts in less than 5 samples
dds_BINF_RvNR_W0_TOF <- dds[,colnames(dds) %in% rownames(BINF_RvNR_W0_TOF_samples)]
dds_BINF_RvNR_W0_TOF_fsub <- dds_BINF_RvNR_W0_TOF[rowSums(counts(dds_BINF_RvNR_W0_TOF) >= 10) >= 5,]
vst_BINF_RvNR_W0_TOF_fsub <- vst(dds_BINF_RvNR_W0_TOF_fsub)

svd_vst_BINF_RvNR_W0_TOF_fsub <- svd(t(assay(vst_BINF_RvNR_W0_TOF_fsub)))
svd_vst_BINF_RvNR_W0_TOF_fsub_df <- data.frame(PC1 = svd_vst_BINF_RvNR_W0_TOF_fsub$u[,1],
                                               PC2 = svd_vst_BINF_RvNR_W0_TOF_fsub$u[,2],
                                               Response = colData(vst_BINF_RvNR_W0_TOF_fsub)$Response,
                                               Week = colData(vst_BINF_RvNR_W0_TOF_fsub)$Week,
                                               DonorID = colData(vst_BINF_RvNR_W0_TOF_fsub)$DonorID,
                                               Sex = colData(vst_BINF_RvNR_W0_TOF_fsub)$Sex,
                                               row.names = colnames(vst_BINF_RvNR_W0_TOF_fsub))

svd_vst_BINF_RvNR_W0_TOF_fsub_var <- svd_vst_BINF_RvNR_W0_TOF_fsub$d^2/(ncol(vst_BINF_RvNR_W0_TOF_fsub)-1)
svd_vst_BINF_RvNR_W0_TOF_fsub_var_explained <- (svd_vst_BINF_RvNR_W0_TOF_fsub_var/sum(svd_vst_BINF_RvNR_W0_TOF_fsub_var))*100

pca_BINF_RvNR_W0_ggplotobj <- ggplot(svd_vst_BINF_RvNR_W0_TOF_fsub_df, aes(x = PC1, y = PC2)) +
  stat_ellipse(alpha = 0.5, linetype = "dashed", aes(col = Response)) +
  geom_point(shape = 21, size = 3, aes(fill = Response)) +
  scale_color_manual(values = response_colors) +
  scale_fill_manual(values = response_colors) +
  labs(# x = paste0("PC1 (", round(svd_vst_BINF_RvNR_W0_TOF_fsub_var_explained[1], 2), "%)"),
       # y = paste0("PC2 (", round(svd_vst_BINF_RvNR_W0_TOF_fsub_var_explained[2], 1), "%)"),
       ) +
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

volcanoplot_degs_BINF_RvNR_W0_TOF_ggplotobj <- degs_BINF_RvNR_W0_TOF %>%
  dplyr::mutate(mlog10padj = -log10(padj),
                Significance = factor(case_when(
                  log2FoldChange > 0 & padj < 0.05 ~ "R",
                  log2FoldChange < 0 & padj < 0.05 ~ "NR",
                  .default = "NS"), 
                  levels = c("R", "NR", "NS"))) %>%
  ggplot(aes(y = mlog10padj, x = log2FoldChange, color = Significance)) +
  geom_hline(yintercept = 0, color = "gray40") +
  geom_vline(xintercept = 0, color = "gray40") +
  geom_point_rast(size = 1.4, alpha = 0.8) +
  geom_label_repel(
    data = (. %>%
              dplyr::filter(gene_symbol %in% goi_figure2$Gene)),
    aes(label = gene_symbol, fill = Significance), 
    nudge_y = 0.5,
    size = 2.8,
    max.overlaps = Inf,
    force = 3, 
    force_pull = 0.5, 
    max.iter = 10000, 
    segment.colour = "black",
    segment.size = 0.6,
    segment.alpha = 0.9,
    point.padding = unit(0.6, "lines"),
    box.padding = unit(0.4, "lines"),
    min.segment.length = 0,
    color = "#000000", # White text
    fontface = "bold", # Bold for contrast
    label.padding = unit(0.15, "lines"),
    show.legend = FALSE) +
  scale_fill_manual(values = response_colors,
                     drop = FALSE) +
  scale_color_manual(values = response_colors,
                     drop = FALSE) +
  labs(x = "log2(fold-change)", 
       y = "-log10(padj)") +
  guides(fill=guide_legend(title="R")) +
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        # plot.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# Cv1

exprs_BINF_RvNR_W0_TOF <- assay(vst_BINF_RvNR_W0_TOF_fsub)
rownames(exprs_BINF_RvNR_W0_TOF) <- gsub("^(ENSG[0-9]+)\\.[0-9]+$", "\\1", rownames(exprs_BINF_RvNR_W0_TOF))
exprs_goi_BINF_RvNR_W0_TOF <- exprs_BINF_RvNR_W0_TOF[goi_figure2$ENSG,]
exprs_scaled_goi_BINF_RvNR_W0_TOF <- t(apply(as.matrix(exprs_goi_BINF_RvNR_W0_TOF), MARGIN = 1, FUN = scale))
colnames(exprs_scaled_goi_BINF_RvNR_W0_TOF) <- colnames(exprs_goi_BINF_RvNR_W0_TOF)

sample_metadata_BINF_RvNR_W0_TOF <- colData(vst_BINF_RvNR_W0_TOF_fsub) %>%
  data.frame()

feature_metadata_goi_BINF_RvNR_W0_TOF <- data.frame(Category = goi_figure2$Category,
                                                    row.names = goi_figure2$ENSG)

exprs_goi_BINF_R_W8vW0_TOF_sample_annotation <- HeatmapAnnotation(
  Response = sample_metadata_BINF_RvNR_W0_TOF$Response,
  col = list(Response = response_colors),
  gp = gpar(col = "#D3D3D3")
)

Heatmap(exprs_scaled_goi_BINF_RvNR_W0_TOF[rownames(feature_metadata_goi_BINF_RvNR_W0_TOF),], 
        column_split = colData(vst_BINF_RvNR_W0_TOF_fsub)$Response,
        top_annotation = exprs_goi_BINF_R_W8vW0_TOF_sample_annotation,
        rect_gp = gpar(col = "#D3D3D3", lwd = 1),
        row_labels = goi_figure2$Gene,
        name = "Expression",
        cluster_columns = T,
        show_row_names = T,
        show_column_names = F,
        row_title = NULL,
        column_title = NULL)

# Cv2
boxplot_BINF_RvNR_W0_TOF_ggplotobj <- data.frame(exprs_BINF_RvNR_W0_TOF[goi_figure2$ENSG,],
                                                 ENSG = goi_figure2$ENSG,
                                                 Gene = goi_figure2$Gene) %>%
  tidyr::pivot_longer(-c(ENSG, Gene), names_to = "Sample_ID", values_to = "Exprs") %>%
  dplyr::left_join(colData(vst_BINF_RvNR_W0_TOF_fsub) %>%
                     data.frame() %>%
                     tibble::rownames_to_column(var = "Sample_ID") %>%
                     dplyr::select(Sample_ID, Response),
                   by = "Sample_ID") %>%
  dplyr::left_join(degs_BINF_RvNR_W0_TOF %>%
                     dplyr::mutate(padj_label = ifelse(padj<0.01, formatC(padj, format = "e", digits = 3), round(padj, 3))),
                   by = c("Gene" = "gene_symbol")) %>%
  dplyr::mutate(Label = paste0(Gene, "\np = ", padj_label),
                Gene = factor(Gene, levels = c("STING1", "IL36G", "IL36RN", "PRKCQ", "NFATC2"))) %>%
  ggplot(aes(x = Response, y = Exprs)) +
  geom_boxplot(aes(fill = Response), alpha = 0.25, outlier.shape = NA) +
  geom_jitter(aes(fill = Response), show.legend = F, shape = 21, width = 0.1) +
  labs(x = "Response",
       y = "Expression") +
  scale_color_manual(values = response_colors) + 
  scale_fill_manual(values = response_colors) + 
  facet_wrap(~Label, scales = "free", nrow = 1) +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

# D

fgsea_BINF_RvNR_W0_TOF <- fgsea_BINF_RvNR_W0_TOF %>%
  dplyr::mutate(Direction = ifelse(NES<0, "NR", "R"),
                Significance = ifelse(padj<0.05, "Significant", "NS"))

fgsea_BINF_RvNR_W0_ggplotobj <- fgsea_BINF_RvNR_W0_TOF %>%
  dplyr::inner_join(pwoi_figure2, by = c("pathway" = "Geneset")) %>%
  dplyr::mutate(aNES = abs(NES)) %>%
  ggplot(aes(x = aNES, y = forcats::fct_reorder(Label, abs(NES)))) +
  geom_bar(aes(fill = Direction, alpha = Figure2, col = Figure2), stat = "identity") +
  scale_alpha_manual(values = c("Highlight" = 1, "Include" = 0.25)) +
  scale_fill_manual(values = response_colors) +
  scale_color_manual(values = c("Include" = "#FFFFFF", "Highlight" = "#000000")) +
  labs(title = "GSEA: Top 15",
       subtitle = "Baseline: R vs NR",
       x = "|NES|") +
  # scale_y_discrete(labels = label_wrap(30)) +
  guides(alpha = "none",
         col = "none") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

# Assembly
figA <- pca_BINF_RvNR_W0_ggplotobj
figB <- volcanoplot_degs_BINF_RvNR_W0_TOF_ggplotobj
figC <- boxplot_BINF_RvNR_W0_TOF_ggplotobj
figD <- fgsea_BINF_RvNR_W0_ggplotobj

fig2 <- figA + 
  figB +
  figC +
  free(figD) +
  plot_layout(design = "
AB
AB
CC
DD
DD
DD
              ")

# Save
ggsave(filename = figA_pdf, 
       plot = figA, 
       width = 5,
       height = 5.5)
ggsave(filename = figB_pdf, 
       plot = figB, 
       width = 5,
       height = 5.5)
ggsave(filename = figC_pdf, 
       plot = figC, 
       width = 7.5,
       height = 2.5)
ggsave(filename = figD_pdf, 
       plot = figD, 
       width = 7.5,
       height = 3.5)