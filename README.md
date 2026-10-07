# ELNID001000000000496617

## Multi-omic analysis of tofacitinib response in ulcerative colitis

Reproducible Snakemake workflow accompanying the manuscript:

> **_[Manuscript title]_**
> [Author list]. _[Journal]_, [year]. DOI: `[TODO: add on acceptance]`

<!-- Internal project identifier: ELNID001000000000496617 -->

This repository contains the complete code (Snakemake workflow, R/Python scripts, conda environments and versions and configuration) used to perform the analyses to characterize the molecular response to tofacitinib in ulcerative colitis (UC), and to compare it against infliximab. The pipeline integrates four data modalities across two tissue compartments, three different studies, and two treatments, harmonizing them into a single reproducible workflow.

The code is released so that the analyses reported in the manuscript can be inspected and, where the underlying data are available, re-run. It is **not** a general-purpose scRNA-seq/RNA-seq/microarray analysis tool. Paths, sample metadata and comparison definitions are specific to this study, but the structure is transparent enough to adapt.

---

## Table of contents

- [Overview](#Overview)
- [Pipeline description](#Pipeline-description)
- [File structure](#file-structure)
- [Software requirements](#software-requirements)
- [Input data](#input-data)
- [Configuration](#configuration)
- [Comparison naming convention](#comparison-naming-convention)
- [Running the workflow](#running-the-workflow)
- [Outputs](#outputs)
- [Notes, limitations and caveats](#notes-limitations-and-caveats)
- [Generative AI usage statement](#generative-ai-usage-statement)
- [Citation](#citation)
- [License and contact](#license-and-contact)

---

## Overview

Inflamed and non-inflamed colonic biopsies and peripheral blood (leukocytes and serum) were collected from UC patients starting tofacitinib treatment. Samples were obtained before treatment (week 0) and at response assessment (week 8). Patients were classified by clinical response (decrease endoscopic Mayo score >= 2 AND decrease Robart Histology Index >= [TODO: FIX]). The transcriptome (bulk RNA-seq) and three targeted proteomic 96-panels (Olink) were profiled, and the findings were placed in context using two external, publicly available cohorts:

| Modality | Treatment | Source | Platform |
|---|---|---|---|
| Bulk RNA-seq | Tofacitinib | This study | RNA-seq, STAR + featureCounts + DESeq2 |
| Proteomics | Tofacitinib | This study | Olink Explore, limma |
| Microarray | Infliximab | Toedter et al. 2011 (GSE23597) | Affymetrix Human Genome U133 Plus 2.0 Array, limma |
| Single-cell RNA-seq (pseudobulk) | Tofacitinib | Melón-Ardanaz et al. 2025 (GSE253006) | 10x Genomics, Seurat pseudobulk + DESeq2 |

Both external datasets differ from the in-house cohort in the clinical metadata available and in their exact definitions of response; these differences are described in the manuscript and the header comments of the corresponding rule files (`workflow/rules/microarray_toedter2011.smk` and `workflow/rules/scrnaseq_melonardanaz2025.smk`) and should be borne in mind when interpreting cross-cohort comparisons. Agreement across these cohorts is treated as cross-platform corroboration.

---

## Pipeline description

Each modality follows the same conceptual path: prepare a data object, subset per comparison, run differential analysis, run gene-set enrichment:

1. **Bulk RNA-seq (own, tofacitinib).** Copy FASTQ and the STAR index to scratch → FastQC → STAR paired-end alignment → filter (remove unmapped / multi-mapped / low-MAPQ / mitochondrial reads) and index with SAMtools → featureCounts → MultiQC → build a `DESeqDataSet` → per-comparison subsetting → DESeq2 differential expression → fgsea.
2. **Olink proteomics (own, tofacitinib).** Build a `SummarizedExperiment` from NPX values → per-comparison subsetting → limma differential abundance → fgsea.
3. **Microarray (Toedter 2011, infliximab).** Import GSE23597 via GEOquery into a `SummarizedExperiment` → subset → limma → fgsea.
4. **scRNA-seq pseudobulk (Melón-Ardanaz 2025, tofacitinib).** Pseudobulk each sample with Seurat → merge → build a `DESeqDataSet` → subset → DESeq2 → fgsea.
5. **Figures.** Figures 1–5 are assembled from the differential-expression and enrichment outputs above.

MSigDB gene sets (via `msigdbr`) are prepared once and shared by every fgsea step.

---

## File structure

The two archives distributed with this repository unpack into the two directories a standard Snakemake project expects. Arrange them like this and run Snakemake **from the project root**:

```
<project root>/
├── workflow/                     # from workflow.zip
│   ├── Snakefile                 # entry point; targets defined in `rule all`
│   ├── rules/                    # one .smk module per modality
│   │   ├── general.smk           #   MSigDB gene-set preparation (shared)
│   │   ├── rnaseq.smk            #   bulk RNA-seq (tofacitinib)
│   │   ├── olink.smk             #   Olink proteomics (tofacitinib)
│   │   ├── microarray_toedter2011.smk    # infliximab microarray
│   │   ├── scrnaseq_melonardanaz2025.smk # scRNA-seq pseudobulk
│   │   ├── figures.smk           #   figures 1–5
│   │   └── machine_learning.smk  #   ML aggregation only (see notes below)
│   ├── scripts/                  # R / R Markdown analysis scripts
│   └── envs/                     # per-rule conda environment definitions
│
├── config/                       # from config.zip
│   ├── config.yaml               # all paths, parameters and comparisons
│   ├── metadata/                 # sample / file manifests (xlsx)
│   └── features_of_interest/     # curated gene and gene-set lists (xlsx)
│
├── resources/                    # external raw data — you provide these
│   ├── toedter2011/              #   GSE23597 series matrix
│   └── melonardanaz2025/         #   GSE253006 matrices, one dir per sample
│
└── output/                       # created by the workflow
```

---

## Software requirements

- **Snakemake** (workflow manager) with the **conda** integration enabled (`--software-deployment-method`). Each rule pins its own environment, so you do not need to  install the scientific stack yourself. Snakemake creates the environments on first run.
- **conda / mamba** to solve those environments.
- A **SLURM** cluster is assumed for the sequencing-heavy steps (the study was run on a BeeGFS-backed HPC running RHEL). The workflow runs locally too, but STAR alignment and pseudobulk merging have substantial memory requirements (see the `resources:` directives in the rules).

Key tool versions (pinned in `workflow/envs/`):

| Tool | Version |
|---|---|
| STAR | 2.7.11b |
| Subread / featureCounts | 2.0.6 |
| SAMtools | 1.19 |
| FastQC / MultiQC | 0.12.1 / 1.21 |
| R Bioconductor DESeq2 | 1.50.2 |
| R Bioconductor limma | 3.66.0 |
| R Bioconductor fgsea | 1.36.2 |
| R Bioconductor SummarizedExperiment | 1.40.0 |
| R Bioconductor GEOquery | 2.78.0 |
| R Seurat | 5.5.1 |
| R msigdbr | 26.1.0 |
| R tidyverse | 2.0.0 |

---

## Input data

The inputs fall into three categories, and this matters for reproducibility:

**1. Provided in this repository** — sample/file manifests and curated
gene/gene-set lists under `config/metadata/` and `config/features_of_interest/`.

**2. Public, must be downloaded into `resources/`:**

- **GSE23597** (Toedter et al. 2011, infliximab microarray). Place the series
  matrix at the path given by `data.toedter2011_gse` in `config.yaml`
  (`resources/toedter2011/GSE23597_series_matrix.txt`).
- **GSE253006** (Melón-Ardanaz et al. 2025, scRNA-seq). Place the per-sample
  count matrices under `data.melonardanaz2025_mtx_dir`
  (`resources/melonardanaz2025/<Sample_ID>/`).

**3. Controlled-access, not included** — the in-house raw sequencing reads
(FASTQ) and Olink NPX measurements are individual-level human data and are not
distributed here. They are available under the terms described in the manuscript's
data-availability statement `[TODO: add repository / accession / access
procedure]`. Without them, the RNA-seq alignment and Olink preparation steps
cannot be re-run, but every downstream differential-expression and enrichment
step can be reproduced from the intermediate objects if those are made available.

**Reference genome/annotation.** The STAR index (Ensembl GRCh38) and GTF
(GENCODE v38) are referenced in `config.yaml` by absolute cluster paths under
`reference:`. Point these at your own local copies before running the alignment
rules.

---

## Configuration

Everything the workflow needs is declared in `config/config.yaml`:

- `base_dir`, `scratch_dir` — project root and fast scratch location on the HPC. Edit these for your environment.
- `R:` — shared helper scripts and the curated gene / gene-set spreadsheets used by the figure scripts.
- `data:` — paths to every metadata manifest and external dataset.
- `reference:` — STAR index and GTF paths.
- `analysis_params:` — filtering and enrichment thresholds (DESeq2 count filter, limma detection filter, fgsea min/max set size and permutations).
- `interaction:` — the a priori gene set used by the tofacitinib-vs-infliximab interaction test.
- Three comparison blocks (see below).

### Comparison blocks

Comparisons are defined declaratively so the same rules apply to every contrast:

- `comparisons` — the in-house tofacitinib cohort (RNA-seq and/or Olink).
- `comparisons_toedter2011` — the infliximab microarray cohort.
- `comparisons_melonardanaz2025` — the tofacitinib scRNA-seq pseudobulk cohort.

Each comparison entry uses these keys:

| Key | Meaning |
|---|---|
| `assays` | Which modalities the comparison applies to (e.g. `[rnaseq, olink]`). Comparisons are routed to the RNA-seq and/or Olink rules accordingly. |
| `filter` | Sample-metadata columns and allowed values used to select the samples for this contrast. |
| `variable` | The tested factor (e.g. `Week` for longitudinal, `Response` for cross-sectional). |
| `contrast` | `[numerator, denominator]` for the effect of interest. |
| `block` | _(optional)_ within-subject blocking factor (`DonorID`) for paired designs. |
| `covariates` | _(optional)_ `numeric:` and `categorical:` adjustment covariates. |

---

## Comparison naming convention

Comparison identifiers encode the compartment, the response group and the
contrast:

```
  BINF_R_W8vW0
  │    │  └── contrast: week 8 vs week 0 (longitudinal)
  │    └───── response group
  └────────── compartment
```

**Compartment:** `BINF` = inflamed biopsy. `BNINF` = non-inflamed biopsy. `PB` = peripheral blood.

**Response group:** `R` = responder. `NR` = non-responder. `LR` = late responder.

**Contrast:** `R_W8vW0` = longitudinal (week 8 vs week 0, within responders). `NR_W8vW0` = longitudinal (week 8 vs week 0, within non-responders). `RvNR_W0` = cross-sectional (responders vs non-responders at baseline).

---

## Running the workflow

From the project root, with conda/mamba and Snakemake available:

```bash
# 1. Test pipeline using a dry run without executing anything
snakemake -np

# 2. Visualize the job DAG
snakemake --dag | dot -Tpdf > dag.pdf

# 3. Run locally with 16 cores, letting Snakemake build the conda environments
snakemake --software-deployment-method conda --cores 16
```

On a SLURM cluster, supply a cluster profile (not included here due to site-specificity) and let Snakemake submit jobs. For example, with the SLURM executor:

```bash
snakemake --software-deployment-method conda --executor slurm --jobs 100 --default-resources slurm_partition="<partition>"
```

The default targets are declared in `rule all` in the `Snakefile`. Many intermediate targets are commented out there so that a full run produces the enrichment tables, the merged pseudobulk matrix and the figures. Uncomment the relevant lines (or name a specific output on the command line) to materialize intermediate objects such as BAMs, count matrices or per-comparison `DESeqDataSet`/`SummarizedExperiment` objects.

---

## Outputs

All results are written under `output/`, organized by modality and comparison, for example:

```
output/
├── general/gs_msigdb.Rds
├── rnaseq/
│   ├── counts/counts.txt
│   ├── multiqc/
│   ├── dds/dds.Rds
│   └── analyses/<comparison>/
│       ├── dds/   (subset DESeqDataSet + VST)
│       ├── degs/  (differential-expression table, .csv + .Rds)
│       └── fgsea/ (enrichment table, .csv)
├── olink/analyses/<comparison>/{se,deps,fgsea}/
├── microarray_toedter2011/analyses/<comparison>/{se,degs,fgsea}/
├── scrnaseq_melonardanaz2025/.../{import,merged,dds,degs,fgsea}/
└── figures/figure{1..5}/*.pdf
```

Differential tables carry effect sizes, standard errors and adjusted p-values. Fgsea tables carry normalized enrichment scores, adjusted p-values andpipe-delimited leading-edge gene lists.

---

## Notes, limitations, and caveats

- **`workflow/scripts/functions.R`** is referenced by `config.yaml`
  (`R.functions`) and loaded by several preparation scripts. Confirm it is present before running; if it is missing from your checkout, add it.
- **External-cohort differences.** The infliximab (Toedter 2011) and scRNA-seq (Melón-Ardanaz 2025) datasets lack some clinical covariates present in the in-house cohort and define response differently. Adjustment models for these cohorts are therefore simpler than for the in-house data, and cross-cohort comparisons are interpreted with that asymmetry in mind.
- **Reproducibility scope.** Steps upstream of the differential-expression tables depend on controlled-access raw data (see [Input data](#input-data)). The enrichment and figure steps are fully reproducible from intermediate objects.

---

## Generative AI usage statement 

Parts of this project were developed with the assistance of Anthropic Claude, Opus 4.8. It was used for generating boilerplate code, debugging, and drafting documentation (including this README.md). All code and documentation has been reviewed, further developed, and tested by a human.

---

## Citation

If you use this code, please cite the manuscript:

```
[TODO: full citation on acceptance — authors, title, journal, year, DOI]
```

and the external datasets it builds on:

- Toedter G, et al. Gene expression profiling and response signatures associated
  with differential responses to infliximab treatment in ulcerative colitis.
  _Am J Gastroenterol_ 2011. GSE23597.
- Melón-Ardanaz E, et al. 2025. _J Crohns Colitis_. DOI: 10.1093/ecco-jcc/jjaf076.
  GSE253006.

---

## License and contact

- **License:** `[TODO: add a LICENSE file — e.g. MIT for code]`
- **Corresponding author / maintainer:** `[TODO: name, ORCID, contact]`