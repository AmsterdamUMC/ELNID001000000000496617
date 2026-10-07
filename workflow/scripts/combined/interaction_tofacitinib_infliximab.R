#!/usr/bin/env Rscript
# Interaction analysis: tofacitinib responders vs infliximab responders over time

###############################################################################
# Purpose
#   Formally test whether the longitudinal (week 8 vs baseline) change in responders DIFFERS between tofacitinib and infliximab, for an a priori non-canonical NF-kB gene set.
#
# Rationale / model
#   Each within-drug longitudinal log2 fold-change is the *time effect* for that drug. The difference between the two drugs' time effects is therefore the drug x time interaction:
#     delta_g = LFC_TOF_g - LFC_IFX_g
#   estimated as a fixed-effect contrast with combined variance:
#     SE_delta = sqrt(SE_TOF^2 + SE_IFX^2)
#     z = delta / SE_delta 
#     p = 2 * pnorm(-|z|)
#
#   Caveat: data from tofacitinib was run on RNA-seq, data from infliximab was run on microarray. Raw log2FC scales are not identical across platforms; microarray compresses dynamic range relative to RNA-seq, which can inflate |delta|. Section 5 below repeats the test after quantile-normalizing the two platforms' LFC distributions to a common scale, as a robustness check.
# References:
#  - Borenstein, Hedges, Higgins, Rothstein (2009, Introduction to Meta‐Analysis)
#  - Ritchie et al. (2015, Nucleic Acids Res 43:e47)
#  - Smyth (2004, Stat Appl Genet Mol Biol 3:Article 3) 
###############################################################################

# Paths
infliximab_degs_csv <- snakemake@input[["infliximab_degs_csv"]]
tofacitinib_degs_csv  <- snakemake@input[["tofacitinib_degs_csv"]]
interaction_degs_csv  <- snakemake@output[["interaction_degs_csv"]]
interaction_nfkbnc_csv <- snakemake@output[["interaction_nfkbnc_csv"]]

# Logs
log_file <- file(log_file, open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

## a priori hypothesis-driven gene set (non-canonical NF-kB)
noncanonical_nfkb <- c("NFKB2","RELB","CD70","LTA","CD40","CD40LG","TNFRSF8","TNFRSF12A")

# Load
tofacitinib_degs <- read.csv(tofacitinib_degs_csv)
infliximab_degs  <- read.csv(infliximab_degs_csv)

## Tofacitinib: DESeq2
tofacitinib_degs_df <- tofacitinib_degs %>%
  dplyr::select(gene_symbol, log2FoldChange, lfcSE) %>%
  dplyr::rename(gene = gene_symbol,
                lfc = log2FoldChange,
                se = lfcSE)%>%
  dplyr::filter(!gene %in% c(""),
                !is.na(gene),
                is.finite(lfc),
                is.finite(se))

## Infliximab: limma
infliximab_degs_df <- infliximab_degs %>%
  dplyr::rename(gene = gene_symbol,
                lfc = log2FoldChange) %>%
  dplyr::mutate(se = lfc/stat) %>%  #recover SE from the moderated t  (t = logFC / SE  ->  SE = |logFC / t|)
  dplyr::filter(!gene %in% c(""),
                !is.na(gene),
                is.finite(lfc),
                is.finite(se)) %>% # Collapse multiple probes per gene to ONE, choosing the highest-expressed probe (max AveExpr). 
  dplyr::group_by(gene) %>%
  dplyr::filter(AveExpr == max(AveExpr)) %>%
  dplyr::select(gene, lfc, se)

## Interaction analysis


qnorm_to_ref <- function(x) {
  # map values to the quantiles of a standard normal by rank
  r <- rank(x, ties.method = "average")
  qnorm((r - 0.5) / length(x))
}

tofactinib_infliximab_degs <- tofacitinib_degs_df %>%
  dplyr::inner_join(infliximab_degs_df, by = "gene", suffix = c("_tof", "_ifx")) %>%
  dplyr::mutate(delta = lfc_tof - lfc_ifx,
                # In addition, put both platforms' LFC distributions on a shared scale to remove the microarray dynamic-range compression, then recompute delta. This yields a scale-robust effect (not a second p-value). Use it to confirm the direction and ranking of the interaction are not artefacts of platform scale.
                lfc_tof_qn = qnorm_to_ref(lfc_tof),
                lfc_ifx_qn = qnorm_to_ref(lfc_ifx),
                delta_qn = lfc_tof_qn - lfc_ifx_qn, 
                delta_se = sqrt(se_tof^2 + se_ifx^2),
                z = delta/delta_se,
                pvalue = 2*pnorm(-abs(z)),
                padj_genomewide = p.adjust(pvalue, method = "BH"),) %>%
  dplyr::arrange(pvalue)

### A priori non-canonical NF-kB set
tofactinib_infliximab_degs_nfkbnc <- tofactinib_infliximab_degs_df %>%
  dplyr::filter(gene %in% noncanonical_nfkb) %>%
  dplyr::mutate(padj_hypdriv = p.adjust(pvalue, method = "BH"))

### Self-contained test: is the mean interaction effect across the set != 0
t.test(tofactinib_infliximab_degs_nfkbnc$delta, mu = 0)
t.test(tofactinib_infliximab_degs_nfkbnc$delta_qn, mu = 0)

# Save
write.csv(tofactinib_infliximab_degs, tofactinib_infliximab_degs_csv)
write.csv(tofactinib_infliximab_degs_nfkbnc, tofactinib_infliximab_degs_nfkbnc_csv)
sink(type = "message")
sink(type = "output")
close(log_file)