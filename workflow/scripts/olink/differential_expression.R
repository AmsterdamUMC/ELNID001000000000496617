#!/usr/bin/env Rscript
# Perform the differential expression per comparison

# Paths
se_rds <- snakemake@input[["se_rds"]]
limma_params <- snakemake@params[["limma_params"]]
comparison_params <- snakemake@params[["comparison"]]
marraylm_rds <- snakemake@output[["marraylm_rds"]]
deps_csv <- snakemake@output[["deps_csv"]]
log_file <- snakemake@log[[1]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Libraries
library(tidyverse)
library(limma)
library(SummarizedExperiment)

sessionInfo()

# Import full se object
olink_se <- readRDS(se_rds)

min_detected <- limma_params$min_detected

# Filter proteins detected in less than 2 samples
message("limma parameters: ")
print(limma_params)
selected_proteins <- rowSums(!is.na(assay(olink_se))) >= min_detected
message("Subset to ", length(which(selected_proteins)), " proteins based on low counts")

se_subset <- olink_se[selected_proteins,]

# Update design to make categorical variables factors and numerical variables numeric

message("Comparison parameters: ")
print(comparison_params)
covariates <- comparison_params$covariates
variable <- comparison_params$variable
contrast <- c(comparison_params$variable, comparison_params$contrast)
block_var <- comparison_params$block
message("variable: ", contrast[1], ". Comparing ", contrast[2], " relative to ", contrast[3])
factor_vars <- c(variable, covariates$numeric, covariates$categorical)
valid_vars <- c()
for(var in factor_vars){
  if(var %in% colnames(colData(se_subset))){
    if(var %in% covariates$numeric) {
      colData(se_subset)[[var]] <- as.numeric(colData(se_subset)[[var]])
      message(paste(var, "range:", paste(c(min(colData(se_subset)[[var]]), max(colData(se_subset)[[var]])), collapse="-")))
      valid_vars[var] <- TRUE
    } else {
      colData(se_subset)[[var]] <- factor(as.character(colData(se_subset)[[var]]))
      if (length(levels(colData(se_subset)[[var]])) > 1){
        message(paste(var, "levels:", paste(levels(colData(se_subset)[[var]]), collapse=", ")))
        valid_vars[var] <- TRUE
      } else {
        message(paste(var, "not enough levels"))
        valid_vars[var] <- FALSE
      }
    }
    if(var != variable){
      varrank <- Matrix::rankMatrix(cbind(colData(se_subset)[[var]], colData(se_subset)[[variable]]))
      if(varrank < 2){
        message(" ", paste(var, "is collinear with", variable))
        valid_vars[var] <- FALSE
      }
    }
  } else {
    message(paste(var, "not found in colData"))
    valid_vars[var] <- FALSE
  }
}
factor_vars <- factor_vars[valid_vars]
message("factor_vars: ", paste(factor_vars, collapse=", "))
design_formula <- as.formula(paste("~", paste(c(0, factor_vars), collapse=" + ")))
design_matrix <- model.matrix(design_formula, data = colData(se_subset))

message("Design formula: ", deparse(design_formula))
contrast_formula <- paste0(contrast[1], contrast[2], "-", contrast[1], contrast[3])
contrast_matrix <- makeContrasts(contrasts = contrast_formula, 
                                 levels = design_matrix)
message("Contrast formula: ", contrast_formula)
# Run limma
if(!is.null(block_var)){
  message("block: ",  block_var)
  block <- factor(colData(se_subset)[[block_var]])
  corfit <- duplicateCorrelation(assay(se_subset), design_matrix, block = block)
  fit_marraylm <- lmFit(assay(se_subset), design = design_matrix, block = block, correlation = corfit$consensus)
} else{
  fit_marraylm <- lmFit(assay(se_subset), design = design_matrix)
}
fit_marraylm <- contrasts.fit(fit = fit_marraylm, contrasts = contrast_matrix)
fit_marraylm <- eBayes(fit_marraylm)

deps <- topTable(fit = fit_marraylm, coef = contrast_formula, adjust = "fdr", number = Inf) %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "OlinkID") %>%
  dplyr::rename(pvalue = P.Value,
                padj = adj.P.Val,
                log2FoldChange = logFC,
                stat = t) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::left_join(rowData(se_subset) %>%
                     data.frame() %>%
                     tibble::rownames_to_column(var = "OlinkID"),
                   by = "OlinkID") %>%
  dplyr::mutate(gene_symbol = Assay) %>%
  dplyr::select(gene_symbol,
                AveExpr,
                log2FoldChange,
                stat,
                pvalue,
                padj,
                OlinkID,
                UniProt, 
                Panel, 
                Panel_Lot_Nr)

# Save
saveRDS(object = fit_marraylm, file = marraylm_rds, compress = "gzip")
write.csv(deps, file = deps_csv, row.names = F)

sessionInfo()

sink(type = "message")
sink(type = "output")
close(log_file)