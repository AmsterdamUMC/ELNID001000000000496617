#!/usr/bin/env python
"""Compile one supplementary table (.xlsx, single sheet) from a pipeline results csv.

The column layout per source type follows the supplementary table legends in config/supplementary_tables.yaml.
Rows keep the order of the pipeline output (sorted by p-value).
The sheet is named <prefix>_<comparison> (e.g. degs_BINF_R_W8vW0).
"""
import sys

import pandas as pd

results_csv = snakemake.input["results_csv"]
table_xlsx = snakemake.output["table_xlsx"]
source = snakemake.params["source"]
comparison = snakemake.params["comparison"]
log_file = snakemake.log[0]

sys.stdout = sys.stderr = open(log_file, "w")

DEGS_RNASEQ = ["gene_symbol", "ENSG", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
FGSEA = ["pathway", "pvalue", "padj", "log2err", "ES", "NES", "size", "leadingEdge"]

layouts = {
    # source: (sheet prefix, output columns, renames applied before selection)
    "rnaseq_degs": ("degs", DEGS_RNASEQ, {}),
    "melonardanaz2025_degs": ("degs", DEGS_RNASEQ, {}),
    "rnaseq_fgsea": ("fgsea", FGSEA, {"pval": "pvalue"}),
    "melonardanaz2025_fgsea": ("fgsea", FGSEA, {"pval": "pvalue"}),
    "toedter2011_fgsea": ("fgsea", FGSEA, {"pval": "pvalue"}),
    "olink_deps": ("deps", ["gene_symbol", "AveExpr", "log2FoldChange", "stat", "pvalue", "padj",
                            "OlinkID", "UniProt", "Panel", "Panel_Lot_Nr"], {}),
    "toedter2011_degs": ("degs", ["gene_symbol", "AveExpr", "log2FoldChange", "stat", "pvalue", "padj",
                                  "ProbeID"], {}),
}

if source not in layouts:
    raise ValueError(f"Unknown source '{source}'")
prefix, columns, renames = layouts[source]

print(f"Source: {source} ({results_csv})")
# Identifiers are read as text so that symbols such as SEPT9 or MARCH1 are written to Excel as strings, never as dates or numbers.
TEXT_COLUMNS = ["gene_symbol", "ENSG", "ENSGv", "pathway", "leadingEdge", "OlinkID", "UniProt", "Panel", "Panel_Lot_Nr", "ProbeID"]
header = pd.read_csv(results_csv, nrows=0).columns
results = pd.read_csv(results_csv, dtype={c: str for c in TEXT_COLUMNS if c in header}).rename(columns=renames)

missing = [c for c in columns if c not in results.columns]
if missing:
    raise ValueError(f"Columns missing from {results_csv}: {missing}. Available: {list(results.columns)}")
dropped = [c for c in results.columns if c not in columns]
if dropped:
    print(f"Dropping columns not in the legend: {dropped}")

table = results[columns]

sheet_name = f"{prefix}_{comparison}"[:31]  # Excel limit
with pd.ExcelWriter(table_xlsx, engine="openpyxl") as writer:
    table.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.sheets[sheet_name].freeze_panes = "A2"

n_sig = int((table["padj"] < 0.05).sum())
print(f"Wrote {len(table)} rows ({n_sig} with padj < 0.05) to {table_xlsx}, sheet '{sheet_name}'")
