#!/usr/bin/env Rscript
# Figure 1
# A) Workflow (created using BioRender)
# B) Boxplots Robarts Histology Index and eMayo score
# C) Scatterplot Robarts Histology Index and eMayo score

# Paths
sample_metadata_xlsx <- snakemake@input[["sample_metadata_xlsx"]] #"config/metadata/samples/260729_sample_metadata.xlsx"
plotting_parameters_r <- snakemake@params[["plotting_parameters"]] #"workflow/scripts/plotting_parameters.R"
figB_pdf <- snakemake@output[["figB_pdf"]]
figC_pdf <- snakemake@output[["figC_pdf"]]
log_file <- snakemake@log[[1]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Libraries
library(tidyverse)
library(ggplot2)
library(patchwork)
library(readxl)

sessionInfo()

# Import
source(plotting_parameters_r)
sample_metadata <- readxl::read_excel(sample_metadata_xlsx)

donor_metadata <- sample_metadata %>%
  dplyr::select(DonorID, Response, Week, wk0_eMayo, wk8_eMayo, wk0_RHI, wk8_RHI) %>%
  dplyr::filter(Week %in% "W8") %>%
  dplyr::mutate(wk0_eMayo = as.numeric(wk0_eMayo),
                wk8_eMayo = as.numeric(wk8_eMayo),
                wk0_RHI = as.numeric(wk0_RHI),
                wk8_RHI = as.numeric(wk8_RHI),
  ) %>%
  unique() 

# B

RHI_eMayo_RvLRvNR_df <- donor_metadata %>%
  dplyr::arrange(DonorID, Week) %>%
  tidyr::pivot_longer(-c(DonorID, Response, Week), names_to = "Variable", values_to = "Score") %>%
  dplyr::select(-Week) %>%
  dplyr::mutate(Timepoint = gsub("wk([0-9])_.+$", "W\\1", Variable),
                Variable = gsub("wk[0-9]_(.+)$", "\\1", Variable),
                Response = factor(Response, levels = c("R", "LR", "NR")))

boxplot_RHI_eMayo_RvLRvNR_ggplotobj <- RHI_eMayo_RvLRvNR_df %>%
  ggplot(aes(x = Timepoint, y = Score, fill = Response)) +
  geom_boxplot(alpha = 0.25, outlier.shape = NA) +
  geom_line(aes(group = DonorID), linetype = "dashed", alpha = 0.75) +
  geom_point() +
  facet_grid(Variable~Response, scale = "free") +
  scale_fill_manual(values = response_colors) +
  labs(title = "Clinical parameters") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# C

RHI_eMayo_RvLRvNR_wide_df <- RHI_eMayo_RvLRvNR_df %>%
  tidyr::pivot_wider(id_cols = c("DonorID", "Response", "Timepoint"), names_from = "Variable", values_from = "Score") %>%
  dplyr::group_by(DonorID) %>%
  dplyr::mutate(deMayo = eMayo[Timepoint == "W8"] - eMayo[Timepoint == "W0"],
                dRHI = RHI[Timepoint == "W8"] - RHI[Timepoint == "W0"]) %>%
  dplyr::select(DonorID, Response, deMayo, dRHI) %>%
  unique()

scatterplot_dRHI_deMayo_RvLRvNR_ggplotobj <- RHI_eMayo_RvLRvNR_wide_df %>%
  ggplot(aes(x = deMayo, y = dRHI, fill = Response)) +
  geom_hline(yintercept = -3, col = "#000000", linetype = "dashed") +
  geom_vline(xintercept = -1, col = "#000000", linetype = "dashed") +
  geom_point(size = 3, shape = 21) +
  scale_fill_manual(values = response_colors) +
  labs(title = "W8-W0") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(face = "bold"))

# Assembly
figB <- boxplot_RHI_eMayo_RvLRvNR_ggplotobj
figC <- scatterplot_dRHI_deMayo_RvLRvNR_ggplotobj

figBC <- figB +
  figC + 
  plot_layout(guides = "collect",
              widths = c(1.5, 1))

# Save

ggsave(filename = figB_pdf, 
       plot = figB, 
       width = 5,
       height = 5)

ggsave(filename = figC_pdf, 
       plot = figC, 
       width = 5.5,
       height = 5)

sink(type = "message")
sink(type = "output")
close(log_file)