#!/usr/bin/env Rscript
# Perform the differential expression per comparison

# Libraries
library(tidyverse)
library(DESeq2)
library(Matrix)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Paths
dds_rds <- snakemake@input[["dds_rds"]]
deseq2_params <- snakemake@params[["deseq2_params"]]
comparison_params <- snakemake@params[["comparison"]]
deseq_rds <- snakemake@output[["deseq_rds"]]
degs_csv <- snakemake@output[["degs_csv"]]
log_file <- snakemake@log[[1]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

# Import full dds object
dds <- readRDS(dds_rds)

min_counts <- deseq2_params$min_counts
min_samples <- deseq2_params$min_samples

# Filter genes with at least  `min_counts` counts in less than `min_samples` samples
message("DESeq2 parameters: ")
print(deseq2_params)
selected_genes <- rowSums(counts(dds) >= min_counts) >= min_samples
message("Subset to ", length(which(selected_genes)), " genes based on filters")

dds_subset <- dds[selected_genes,]

# Update design to make categorical variables factors and numerical variables numeric

message("Comparison parameters: ")
print(comparison_params)
covariates <- comparison_params$covariates
variable <- comparison_params$variable
contrast <- c(comparison_params$variable, comparison_params$contrast)
message("variable: ", contrast[1], ". Comparing ", contrast[2], " relative to ", contrast[3])
factor_vars <- c(variable, covariates$numeric, covariates$categorical)
valid_vars <- c()
for(var in factor_vars){
  if(var %in% colnames(colData(dds_subset))){
    if(var %in% covariates$numeric) {
      colData(dds_subset)[[var]] <- as.numeric(colData(dds_subset)[[var]])
      message(paste(var, "range:", paste(c(min(colData(dds_subset)[[var]]), max(colData(dds_subset)[[var]])), collapse="-")))
      valid_vars[var] <- TRUE
    } else {
      colData(dds_subset)[[var]] <- factor(as.character(colData(dds_subset)[[var]]))
      if (length(levels(colData(dds_subset)[[var]])) > 1){
        message(paste(var, "levels:", paste(levels(colData(dds_subset)[[var]]), collapse=", ")))
        valid_vars[var] <- TRUE
      } else {
        message(paste(var, "not enough levels"))
        valid_vars[var] <- FALSE
      }
    }
    if(var != variable){
      varrank <- Matrix::rankMatrix(cbind(colData(dds_subset)[[var]], colData(dds_subset)[[variable]]))
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
design_formula <- as.formula(paste("~", paste(factor_vars, collapse=" + ")))

design(dds_subset) <- design_formula

message("Design formula: ", deparse(design_formula))

# Run DESeq2
dds_subset <- DESeq(dds_subset)
print(results(dds_subset, contrast = contrast))
degs <- results(dds_subset, contrast = contrast) %>%
  data.frame(., gene_symbol = rownames(.)) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::mutate(ENSG = mapIds(org.Hs.eg.db,
                              keys = gene_symbol,
                              column = "ENSEMBL",
                              keytype = "SYMBOL",
                              multiVals = "first")) %>%
  dplyr::select(gene_symbol,
                ENSG,
                baseMean,
                log2FoldChange,
                lfcSE,
                stat,
                pvalue,
                padj)

# Save
saveRDS(object = dds_subset, file = deseq_rds, compress = "gzip")
write.csv(degs, file = degs_csv, row.names = F)

sessionInfo()

sink(type = "message")
sink(type = "output")
close(log_file)