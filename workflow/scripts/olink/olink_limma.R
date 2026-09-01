#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(OlinkAnalyze)
  library(limma)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ggrepel)
})

# ---------- Open the log file for recording ----------
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, append = FALSE, type = "output")   # capture stdout
sink(log_file, append = TRUE,  type = "message")  # capture stderr / message()

# ---------- Access snakemake inputs and params ----------
message("--- Loading Snakemake ---")
source("workflow/scripts/functions.R")
mtx_file       <- snakemake@input[["mtx_file"]]
meta_file      <- snakemake@input[["meta_file"]]
feature_file   <- snakemake@input[["feature_file"]]
filter         <- snakemake@params$filter
limma_opts     <- snakemake@params$limma_opts
output_csv     <- snakemake@output[["results_csv"]]

# ---------- Loading data ----------
message("--- Loading OLINK dataset ---")
olink_mtx <- read.csv(mtx_file, row.names=1)
meta <- readxl::read_excel(meta_file)
meta <- as.data.frame(meta)
feature_meta <- read.csv(feature_file, row.names=1)
rownames(meta) <- meta$ClientAccessionID
rownames(meta) <- gsub("-", ".", rownames(meta))
rownames(meta) <- paste0("X", rownames(meta))
message("Original samples: ", ncol(olink_mtx), "\tOriginal OLINK features: ", nrow(olink_mtx))

# ---------- Subset samples ----------
message("--- Subset Samples ---")
# Filter based on filters
keep <- rep(TRUE, nrow(meta))
for (col in names(filter)) {
  if(col!="Status"){
    vals <- filter[[col]]
    message(paste(c(col, paste(vals)), collapse = "->"))
    keep <- keep & meta[[col]] %in% vals
  }
}
present_samples <- intersect(names(olink_mtx), rownames(meta[keep,]))
olink_mtx <- olink_mtx[,present_samples]
meta <- meta[present_samples,]
message("Subset to ", ncol(olink_mtx), " samples based on filters")
# Filter proteins detected in less than 2 samples
keep <- rowSums(!is.na(olink_mtx)) >= 2
olink_mtx <- olink_mtx[keep, ]
message("Subset to ", nrow(olink_mtx), " OLINK features based on low counts")

# ---------- Update design ----------
message("--- limma Experiment Deisgn ---")
confounders    <- limma_opts$design$confounders
variable      <- limma_opts$design$variable
block       <- limma_opts$design$block
factor_vars <- paste(c(confounders$numeric, confounders$categorical, variable))
valid_vars <- c()
for (var in factor_vars) {
  if (var %in% colnames(meta)) {
    if (var %in% confounders$numeric) {
      meta[[var]] <- as.numeric(meta[[var]])
      message(paste(var, "range:", paste(c(min(meta[[var]]), max(meta[[var]])), collapse="-")))
      valid_vars[var] <- TRUE
    } else {
      meta[[var]] <- factor(clean_string(as.character(meta[[var]])))
      if (length(levels(meta[[var]])) > 1){
        message(paste(var, "levels:", paste(levels(meta[[var]]), collapse=", ")))
        valid_vars[var] <- TRUE
      } else {
        message(paste(var, "not enough levels"))
        valid_vars[var] <- FALSE
      }
    }
  } else {
    message(paste(var, "not found in colData"))
    valid_vars[var] <- FALSE
  }
}
factor_vars <- factor_vars[valid_vars]
design_formula <- as.formula(
  paste("~", paste(factor_vars, collapse=" + "))
)
design <- model.matrix(design_formula, data = meta)
message(paste("Design formula: ", design_formula))

# ---------- Run limma ----------
if (!is.null(block)){
    message("block: ",  block)
    block <- factor(meta[[block]])
    # ---------- Run limma ----------
    corfit <- duplicateCorrelation(olink_mtx, design, block = block)
    fit <- lmFit(olink_mtx, design, block = block, correlation = corfit$consensus)
    fit <- eBayes(fit)
} else{
    fit <- lmFit(olink_mtx, design)
    fit <- eBayes(fit)
}
#dep <- topTable(fit, coef=paste0(variable, levels(meta[[variable]])[2]),adjust="fdr", number=Inf)
if (all(limma_opts$contrast == levels(meta[[variable]]))){
    contrast_formula <- paste0("-", paste0(variable, levels(meta[[variable]])[2])) 
} else if(all(rev(limma_opts$contrast) == levels(meta[[variable]]))){
    contrast_formula <- paste0(variable, levels(meta[[variable]])[2])
}
contrast_matrix <- makeContrasts(
  contrasts = contrast_formula,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
dep <- topTable(fit2, adjust.method = "fdr", number = Inf)
dep$gene_symbol <- feature_meta[row.names(dep),"Assay"]
colnames(dep) <- c('log2FoldChange', 'AveExpr',	't', 'pvalue',	'padj', 'B', 'gene_symbol')

dir.create(dirname(output_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(dep, output_csv, row.names = FALSE)
message("limma results saved: ", output_csv)

sink(type = "message")
sink(type = "output")
close(log_file)



