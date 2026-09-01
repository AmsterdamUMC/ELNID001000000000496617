###############################################################################
# Drug x time interaction test: tofacitinib vs infliximab (responders)
# -----------------------------------------------------------------------------
# Purpose
#   Formally test whether the longitudinal (week 8 vs baseline) change in
#   responders DIFFERS between tofacitinib and infliximab, for an a priori
#   non-canonical NF-kB gene set.
#
# Rationale / model
#   Each within-drug longitudinal log2 fold-change is the *time effect* for
#   that drug. The difference between the two drugs' time effects is therefore
#   the drug x time interaction:
#         delta_g = LFC_tofa,g  -  LFC_ifx,g
#   estimated as a fixed-effect contrast with combined variance:
#         SE_delta = sqrt(SE_tofa^2 + SE_ifx^2)
#         z        = delta / SE_delta ,   p = 2 * pnorm(-|z|)
#   A joint per-sample model with an explicit drug:time term would be ideal,
#   but the two arms are on different platforms (RNA-seq/DESeq2 vs Illumina
#   microarray/limma) and only summary statistics are available, so they
#   cannot be co-modelled at the observation level. This contrast is the
#   correct summary-statistic equivalent of that interaction.
#
# CAVEAT (report this): raw log2FC scales are not identical across platforms;
#   microarray compresses dynamic range relative to RNA-seq, which can inflate
#   |delta|. Section 5 below repeats the test after quantile-normalising the
#   two platforms' LFC distributions to a common scale, as a robustness check.
###############################################################################

## ------------------------------------------------------------------ 0. paths
tofa_path <- "rna_Tissue_R_Inflamed_Week_8_vs_Baseline_deseq2_results.csv"
ifx_path  <- "Infliximab_R_Week_8_vs_Baseline_limma_results.csv"
out_path  <- "drug_time_interaction_results.csv"

## a priori hypothesis-driven gene set (non-canonical NF-kB)
noncanonical_nfkb <- c("NFKB2","RELB","CD70","LTA","CD40","CD40LG",
                       "TNFRSF8","TNFRSF12A")

## --------------------------------------------------------- 1. load + tidy
tofa <- read.csv(tofa_path, stringsAsFactors = FALSE)
ifx  <- read.csv(ifx_path,  stringsAsFactors = FALSE)

# DESeq2 already provides the standard error of the LFC (lfcSE)
tofa <- data.frame(
  gene   = tofa$gene_symbol,
  lfc_t  = as.numeric(tofa$log2FoldChange),
  se_t   = as.numeric(tofa$lfcSE),
  stringsAsFactors = FALSE
)
tofa <- tofa[!is.na(tofa$gene) & tofa$gene != "" &
               is.finite(tofa$lfc_t) & is.finite(tofa$se_t), ]

# limma: recover SE from the moderated t  (t = logFC / SE  ->  SE = |logFC / t|)
ifx$lfc_i <- as.numeric(ifx$log2FoldChange)
ifx$t     <- as.numeric(ifx$t)
ifx$AveExpr <- as.numeric(ifx$AveExpr)
ifx$se_i  <- abs(ifx$lfc_i / ifx$t)
ifx <- ifx[is.finite(ifx$lfc_i) & is.finite(ifx$se_i) & is.finite(ifx$AveExpr), ]

# collapse multiple probes per gene to ONE, choosing the highest-expressed
# probe (max AveExpr). This is unbiased with respect to the contrast, unlike
# picking the most significant probe. Alternative choices can be swapped here.
ifx <- ifx[order(-ifx$AveExpr), ]
ifx <- ifx[!duplicated(ifx$gene_symbol), ]
ifx <- data.frame(
  gene  = ifx$gene_symbol,
  lfc_i = ifx$lfc_i,
  se_i  = ifx$se_i,
  stringsAsFactors = FALSE
)

## ------------------------------------------------- 2. interaction contrast
m <- merge(tofa, ifx, by = "gene")
m$delta    <- m$lfc_t - m$lfc_i                 # drug x time interaction effect
m$se_delta <- sqrt(m$se_t^2 + m$se_i^2)
m$z        <- m$delta / m$se_delta
m$p        <- 2 * pnorm(-abs(m$z))
m$padj_genomewide <- p.adjust(m$p, method = "BH")   # stringent, all shared genes

cat(sprintf("Genes in interaction model: %d\n\n", nrow(m)))

## --------------------------------------- 3. a priori non-canonical NF-kB set
set <- m[m$gene %in% noncanonical_nfkb, ]
set$padj_withinset <- p.adjust(set$p, method = "BH")  # FDR within the a priori set
set <- set[order(set$delta), ]

cat("Drug x time interaction (Tofa - IFX), non-canonical NF-kB set:\n")
print(format(
  set[, c("gene","lfc_t","lfc_i","delta","z","p",
          "padj_withinset","padj_genomewide")],
  digits = 3), row.names = FALSE)

## ----------------------------------------------------- 4. set-level roll-up
# Self-contained test: is the mean interaction effect across the set != 0?
tt <- t.test(set$delta, mu = 0)
n_sig <- sum(set$padj_genomewide < 0.05 & set$delta < 0)
cat(sprintf(
  "\nSet-level: mean delta = %+.2f  (one-sample t p = %.2e)\n", 
  mean(set$delta), tt$p.value))
cat(sprintf(
  "Genes with tofacitinib significantly MORE downregulated (genome-wide FDR<0.05): %d/%d\n",
  n_sig, nrow(set)))

## --------------------------- 5. sensitivity: common-scale (quantile) check
# Put both platforms' LFC distributions on a shared scale to remove the
# microarray dynamic-range compression, then recompute delta. This yields a
# scale-robust effect (not a second p-value); use it to confirm the DIRECTION
# and ranking of the interaction are not artefacts of platform scale.
qnorm_to_ref <- function(x) {
  # map values to the quantiles of a standard normal by rank (INT)
  r <- rank(x, ties.method = "average")
  qnorm((r - 0.5) / length(x))
}
m$lfc_t_qn <- qnorm_to_ref(m$lfc_t)
m$lfc_i_qn <- qnorm_to_ref(m$lfc_i)
m$delta_qn <- m$lfc_t_qn - m$lfc_i_qn

set_qn <- m[m$gene %in% noncanonical_nfkb, c("gene","delta","delta_qn")]
set_qn <- set_qn[order(set_qn$delta_qn), ]
cat("\nScale-robustness (quantile-normalised delta, same-scale sensitivity):\n")
print(format(set_qn, digits = 3), row.names = FALSE)
cat(sprintf("Sign concordance raw vs quantile-normalised: %d/%d genes\n",
            sum(sign(set_qn$delta) == sign(set_qn$delta_qn)), nrow(set_qn)))

## ------------------------------------------------------------- 6. write out
write.csv(set[, c("gene","lfc_t","se_t","lfc_i","se_i","delta","se_delta",
                  "z","p","padj_withinset","padj_genomewide")],
          out_path, row.names = FALSE)
cat(sprintf("\nWritten: %s\n", out_path))