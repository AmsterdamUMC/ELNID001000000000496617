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

# Pairing / blocking
# If a block (e.g. DonorID) is given AND the tested variable varies within blocks (longitudinal design), the block is added as a fixed effect: ~ block + covariates + variable.
# Covariates that are constant within every block (Age, Sex, Prednisone, Anti_TNF, Vedo, ...) are fully absorbed by the block and would make the model matrix rank-deficient, so they are dropped. Their effect is already controlled for by the within-donor comparison.
# Blocks lacking either contrast level contribute nothing to the contrast and are removed so that only complete pairs are analysed.
# For cross-sectional comparisons (one sample per block) the block is ignored.
block_var <- unlist(comparison_params$block)
if(!is.null(block_var)){
  block_var <- block_var[1]
  if(!block_var %in% colnames(colData(dds_subset))){
    stop("Block variable '", block_var, "' not found in colData")
  }
  block <- factor(as.character(colData(dds_subset)[[block_var]]))
  variable_values <- as.character(colData(dds_subset)[[variable]])
  varies_within_block <- any(tapply(variable_values, block, function(x) length(unique(x)) > 1))

  if(!varies_within_block){
    message("Block '", block_var, "': '", variable, "' does not vary within blocks (cross-sectional); block ignored")
    block_var <- NULL
  } else {
    # Keep complete pairs only
    complete_blocks <- names(which(tapply(variable_values, block, function(x) all(contrast[2:3] %in% x))))
    keep_samples <- as.character(block) %in% complete_blocks
    if(any(!keep_samples)){
      message("Removing ", sum(!keep_samples), " sample(s) from blocks without both '", contrast[2], "' and '", contrast[3], "': ",
              paste(colnames(dds_subset)[!keep_samples], collapse = ", "))
    }
    dds_subset <- dds_subset[, keep_samples]
    colData(dds_subset)[[block_var]] <- droplevels(factor(as.character(colData(dds_subset)[[block_var]])))
    for(var in factor_vars){
      if(is.factor(colData(dds_subset)[[var]])) colData(dds_subset)[[var]] <- droplevels(colData(dds_subset)[[var]])
    }
    message("Block '", block_var, "': ", nlevels(colData(dds_subset)[[block_var]]), " complete pairs, ", ncol(dds_subset), " samples")

    # Drop covariates nested within the block
    block <- colData(dds_subset)[[block_var]]
    candidate_covariates <- setdiff(factor_vars, variable)
    nested <- candidate_covariates[vapply(candidate_covariates, function(var){
      all(tapply(as.character(colData(dds_subset)[[var]]), block, function(x) length(unique(x)) == 1))
    }, logical(1))]
    if(length(nested) > 0){
      message("Dropping covariates constant within '", block_var, "' (absorbed by block): ", paste(nested, collapse = ", "))
    }
    factor_vars <- setdiff(factor_vars, nested)
  }
}

# Tested variable last, preceded by block and remaining covariates
design_vars <- c(block_var, setdiff(factor_vars, variable), variable)
design_formula <- as.formula(paste("~", paste(design_vars, collapse=" + ")))

# Fail early with an informative message rather than inside DESeq()
design_rank <- Matrix::rankMatrix(model.matrix(design_formula, data = as.data.frame(colData(dds_subset))))
design_ncol <- ncol(model.matrix(design_formula, data = as.data.frame(colData(dds_subset))))
if(design_rank < design_ncol){
  stop("Design matrix is not full rank (rank ", design_rank, " < ", design_ncol, " columns): ", deparse(design_formula))
}

design(dds_subset) <- design_formula

message("Design formula: ", deparse(design_formula))

# Run DESeq2
dds_subset <- DESeq(dds_subset)
degs <- results(dds_subset, contrast = contrast) %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "ENSGv") %>%
  dplyr::arrange(pvalue) %>%
  dplyr::mutate(ENSG = gsub("(^ENSG[0-9]+)\\.[0-9]+$", "\\1", ENSGv),
                gene_symbol = mapIds(org.Hs.eg.db,
                                     keys = ENSG,
                                     column = "SYMBOL",
                                     keytype = "ENSEMBL",
                                     multiVals = "first")) %>%
  dplyr::select(gene_symbol,
                ENSG,
                ENSGv,
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