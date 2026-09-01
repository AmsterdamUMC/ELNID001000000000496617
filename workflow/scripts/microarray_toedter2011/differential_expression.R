#!/usr/bin/env Rscript
# Perform the differential expression per comparison

# Paths
se_rds <- snakemake@input[["se_rds"]]
comparison_params <- snakemake@params[["comparison"]]
marraylm_rds <- snakemake@output[["marraylm_rds"]]
degs_csv <- snakemake@output[["degs_csv"]]
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
toedter2011_se <- readRDS(se_rds)

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
  if(var %in% colnames(colData(toedter2011_se))){
    if(var %in% covariates$numeric) {
      colData(toedter2011_se)[[var]] <- as.numeric(colData(toedter2011_se)[[var]])
      message(paste(var, "range:", paste(c(min(colData(toedter2011_se)[[var]]), max(colData(toedter2011_se)[[var]])), collapse="-")))
      valid_vars[var] <- TRUE
    } else {
      colData(toedter2011_se)[[var]] <- factor(as.character(colData(toedter2011_se)[[var]]))
      if (length(levels(colData(toedter2011_se)[[var]])) > 1){
        message(paste(var, "levels:", paste(levels(colData(toedter2011_se)[[var]]), collapse=", ")))
        valid_vars[var] <- TRUE
      } else {
        message(paste(var, "not enough levels"))
        valid_vars[var] <- FALSE
      }
    }
    if(var != variable){
      varrank <- Matrix::rankMatrix(cbind(colData(toedter2011_se)[[var]], colData(toedter2011_se)[[variable]]))
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
design_matrix <- model.matrix(design_formula, data = colData(toedter2011_se))

message("Design formula: ", deparse(design_formula))
contrast_formula <- paste0(contrast[1], contrast[2], "-", contrast[1], contrast[3])
contrast_matrix <- makeContrasts(contrasts = contrast_formula, 
                                 levels = design_matrix)
message("Contrast formula: ", contrast_formula)
# Run limma
if(!is.null(block_var)){
  message("block: ",  block_var)
  block <- factor(colData(toedter2011_se)[[block_var]])
  corfit <- duplicateCorrelation(assay(toedter2011_se, i = "logqnorm"), design_matrix, block = block)
  fit_marraylm <- lmFit(assay(toedter2011_se, i = "log"), design = design_matrix, block = block, correlation = corfit$consensus)
} else{
  fit_marraylm <- lmFit(assay(toedter2011_se, i = "logqnorm"), design = design_matrix)
}
fit_marraylm <- contrasts.fit(fit = fit_marraylm, contrasts = contrast_matrix)
fit_marraylm <- eBayes(fit_marraylm)

degs <- topTable(fit = fit_marraylm, coef = contrast_formula, adjust = "fdr", number = Inf) %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "ProbeID") %>%
  dplyr::rename(pvalue = P.Value,
                padj = adj.P.Val,
                log2FoldChange = logFC,
                stat = t) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::left_join(rowData(toedter2011_se) %>%
                     data.frame() %>%
                     tibble::rownames_to_column(var = "ProbeID"),
                   by = "ProbeID") %>%
  dplyr::mutate(gene_symbol = Gene.Symbol) %>%
  dplyr::select(gene_symbol,
                AveExpr,
                log2FoldChange,
                stat,
                pvalue,
                padj,
                ProbeID)

# Save
saveRDS(object = fit_marraylm, file = marraylm_rds, compress = "gzip")
write.csv(degs, file = degs_csv, row.names = F)

sink(type = "message")
sink(type = "output")
close(log_file)