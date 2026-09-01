matrix

#### Own

repeat_samples_tof <- colData(vst_R_W8vW0_TOF) %>%
  data.frame() %>%
  dplyr::filter(Visit %in% c("Baseline", "Week_8")) %>%
  dplyr::group_by(ClientGroupName.x) %>%
  dplyr::summarize(n = n()) %>%
  dplyr::filter(n>1)

exprs_goi_w8mw0_R_W8vW0_TOF <- assay(vst_R_W8vW0_TOF) %>%
  data.frame() %>%
  dplyr::mutate(ensg = rownames(vst_R_W8vW0_TOF),
                gene_name = mapIds(org.Hs.eg.db, 
                                   keys = ensg,
                                   column = "SYMBOL", 
                                   keytype = "ENSEMBL")) %>%
  dplyr::filter(gene_name %in% goi_figure4$Gene) %>%
  tidyr::pivot_longer(-c(ensg, gene_name), names_to = "Sample_ID", values_to = "Expr") %>%
  dplyr::group_by(ensg, gene_name, Sample_ID) %>%
  dplyr::left_join(colData(vst_R_W8vW0_TOF) %>% 
                     data.frame() %>%
                     dplyr::mutate(ClientAccessionID = paste0("X", ClientAccessionID)) %>%
                     dplyr::select(ClientAccessionID, ClientGroupName.x, Visit, Response),
                   by = c("Sample_ID" = "ClientAccessionID")) %>%
  dplyr::rename(Timepoint = Visit,
                Donor_ID = ClientGroupName.x) %>%
  dplyr::mutate(Timepoint = ifelse(Timepoint == "Baseline", "W0", "W8")) %>%
  dplyr::arrange(gene_name, Donor_ID, Timepoint) %>%
  dplyr::group_by(gene_name, Donor_ID) %>%
  dplyr::mutate(w8mw0 = Expr - Expr[Timepoint == 'W0'])

exprs_goi_w8mw0_R_W8vW0_TOF_wide <- exprs_goi_w8mw0_R_W8vW0_TOF %>%
  dplyr::filter(Timepoint == "W8") %>%
  tidyr::pivot_wider(id_cols = c(gene_name), names_from = Donor_ID, values_from = w8mw0) %>%
  tibble::column_to_rownames(var = "gene_name")

exprs_goi_w8mw0_R_W8vW0_TOF_sample_metadata <- exprs_goi_w8mw0_R_W8vW0_TOF %>%
  dplyr::ungroup() %>%
  dplyr::select(Donor_ID, Response) %>%
  unique()

#### Toedter 2011

repeat_samples_toedter2011 <- pData(eset_RvNR_IFX_toedter2011) %>%
  data.frame() %>%
  dplyr::filter(time.ch1 %in% c("W0", "W8")) %>%
  dplyr::group_by(characteristics_ch1.3) %>%
  dplyr::summarize(n = n()) %>%
  dplyr::filter(n>1)

selected_samples_toedter2011 <- pData(eset_RvNR_IFX_toedter2011) %>%
  data.frame() %>%
  dplyr::filter(dose.ch1 %in% c("5mg/kg", "10mg/kg"),
                time.ch1 %in% c("W0", "W8"),
                wk8.response.ch1 %in% c("Yes"),
                characteristics_ch1.3 %in% repeat_samples_toedter2011$characteristics_ch1.3)

eset_R_W8vW0_IFX_toedter2011 <- eset_RvNR_IFX_toedter2011[,selected_samples_toedter2011$geo_accession]

exprs_goi_w8mw0_R_W8vW0_IFX_toedter2011 <- exprs(eset_R_W8vW0_IFX_toedter2011) %>%
  data.frame() %>%
  dplyr::mutate(gene_name = fData(eset_R_W8vW0_IFX_toedter2011)$`Gene Symbol`) %>%
  dplyr::filter(gene_name %in% goi_figure4$Gene) %>%
  tidyr::pivot_longer(-gene_name, names_to = "Sample_ID", values_to = "Expr") %>%
  dplyr::filter(Sample_ID != "GSM578741") %>% #Replicate
  dplyr::group_by(gene_name, Sample_ID) %>%
  dplyr::summarize(Expr = mean(log10(Expr))) %>%
  dplyr::left_join(pData(eset_R_W8vW0_IFX_toedter2011) %>% 
                     data.frame() %>%
                     dplyr::select(geo_accession, characteristics_ch1.2, characteristics_ch1.3, characteristics_ch1.4),
                   by = c("Sample_ID" = "geo_accession")) %>%
  dplyr::rename(Timepoint = characteristics_ch1.2,
                Donor_ID = characteristics_ch1.3,
                Response = characteristics_ch1.4) %>%
  dplyr::mutate(Timepoint = gsub("time: ", "", Timepoint),
                Donor_ID = gsub("subject: ", "", Donor_ID),
                Response = ifelse(Response == "wk8 response: No", "NR", "R")) %>%
  dplyr::arrange(gene_name, Donor_ID, Timepoint) %>%
  dplyr::group_by(gene_name, Donor_ID) %>%
  dplyr::mutate(w8mw0 = Expr - Expr[Timepoint == 'W0'])

exprs_goi_w8mw0_R_W8vW0_IFX_toedter2011_wide <- exprs_goi_w8mw0_R_W8vW0_IFX_toedter2011 %>%
  dplyr::filter(Timepoint == "W8") %>%
  tidyr::pivot_wider(id_cols = c(gene_name), names_from = Donor_ID, values_from = w8mw0) %>%
  tibble::column_to_rownames(var = "gene_name")

exprs_goi_w8mw0_R_W8vW0_IFX_toedter2011_sample_metadata <- exprs_goi_w8mw0_R_W8vW0_IFX_toedter2011 %>%
  dplyr::ungroup() %>%
  dplyr::select(Donor_ID, Response) %>%
  unique()


exprs_goi_w8mw0_R_W8vW0_TOFvIFX_wide <- data.frame(exprs_goi_w8mw0_R_W8vW0_TOF_wide) %>%
  tibble::rownames_to_column(var = "gene_name") %>%
  dplyr::inner_join(data.frame(exprs_goi_w8mw0_R_W8vW0_IFX_toedter2011_wide) %>%
                      tibble::rownames_to_column(var = "gene_name"),
                    by = "gene_name") %>%
  tibble::column_to_rownames(var = "gene_name")

exprs_goi_w8mw0_R_W8vW0_TOFvIFX_sample_metadata <- exprs_goi_w8mw0_R_W8vW0_TOF_sample_metadata %>%
  dplyr::mutate(Donor_ID = paste0("X", as.character(Donor_ID)),
                Treatment = "Tofacitinib",
                Source = "Own") %>% 
  dplyr::rows_append(exprs_goi_w8mw0_R_W8vW0_IFX_toedter2011_sample_metadata %>%
                       dplyr::mutate(Treatment = "Infliximab",
                                     Source = "Toedter 2011")) %>%
  dplyr::select(Donor_ID, Treatment) %>%
  tibble::column_to_rownames(var = "Donor_ID")

exprs_goi_w8mw0_R_W8vW0_TOFvIFX_feature_metadata <- goi_figure4 %>%
  dplyr::filter(Gene %in% rownames(exprs_goi_w8mw0_R_W8vW0_TOFvIFX_wide)) %>%
  tibble::column_to_rownames(var = "Gene")
exprs_goi_w8mw0_R_W8vW0_TOFvIFX_feature_metadata <- exprs_goi_w8mw0_R_W8vW0_TOFvIFX_feature_metadata[rownames(exprs_goi_w8mw0_R_W8vW0_TOF_wide),]

pheatmap::pheatmap(t(exprs_goi_w8mw0_R_W8vW0_TOFvIFX_wide[,rownames(exprs_goi_w8mw0_R_W8vW0_TOFvIFX_sample_metadata)]),
                   annotation_row = exprs_goi_w8mw0_R_W8vW0_TOFvIFX_sample_metadata, 
                   scale = "column", 
                   cluster_rows = F)

exprs_goi_w8mw0_R_W8vW0_TOFvIFX_gene_annotation <- HeatmapAnnotation(
  Geneset = exprs_goi_w8mw0_R_W8vW0_TOFvIFX_feature_metadata$Label, 
  col = list(Geneset = geneset_label_colors),
  gp = gpar(col = "#D3D3D3")
)

exprs_goi_w8mw0_R_W8vW0_TOFvIFX_sample_annotation <- rowAnnotation(
  Treatment = exprs_goi_w8mw0_R_W8vW0_TOFvIFX_sample_metadata$Treatment, 
  col = list(Geneset = treatment_colors),
  gp = gpar(col = "#D3D3D3")
)

Heatmap(t(exprs_goi_w8mw0_R_W8vW0_TOFvIFX_wide[,rownames(exprs_goi_w8mw0_R_W8vW0_TOFvIFX_sample_metadata)]), 
        name = "log2(W8-W0)", 
        row_split = exprs_goi_w8mw0_R_W8vW0_TOFvIFX_sample_metadata$Treatment,
        column_split = exprs_goi_w8mw0_R_W8vW0_TOFvIFX_feature_metadata$Label,
        top_annotation = exprs_goi_w8mw0_R_W8vW0_TOFvIFX_gene_annotation,
        left_annotation = exprs_goi_w8mw0_R_W8vW0_TOFvIFX_sample_annotation,
        rect_gp = gpar(col = "#D3D3D3", lwd = 1),
        show_row_names = FALSE,
        row_title = NULL,
        column_title = NULL)