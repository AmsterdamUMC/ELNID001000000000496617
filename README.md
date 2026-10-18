# ELNID001000000000496617

## Multi-omic analysis of tofacitinib response in ulcerative colitis

Reproducible Snakemake workflow accompanying the manuscript:

> **_Longitudinal multi-omics profiling of tofacitinib response identifies a pre-existing Janus kinase-tractable inflammatory state in ulcerative colitis_**
> Fung ICN, de Jong DC, Liu K, _et al._ Preprint: medRxiv, DOI: [10.64898/2026.10.07.26364964](https://doi.org/10.64898/2026.10.07.26364964). Journal version: `[TODO: add on acceptance]`

<!-- Internal project identifier: ELNID001000000000496617 -->

This repository contains the complete code (Snakemake workflow, R/Python scripts, conda environments and versions and configuration) used to perform the analyses to characterize the molecular response to tofacitinib in ulcerative colitis (UC), and to compare it against infliximab. The pipeline integrates two data modalities (transcriptomics and proteomics) across mucosal and peripheral compartments, three different studies, and two treatments, harmonizing them into a single reproducible workflow.

The code is released so that the analyses reported in the manuscript can be inspected and, where the underlying data are available, re-run. It is **not** a general-purpose scRNA-seq/RNA-seq/microarray analysis tool. Paths, sample metadata and comparison definitions are specific to this study, but the structure is transparent enough to adapt.

---

## Table of contents

- [Overview](#overview)
- [Pipeline description](#pipeline-description)
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

Inflamed and non-inflamed colonic biopsies, peripheral whole blood and serum were collected from patients with moderate-to-severe UC starting tofacitinib (TOFA-histology trial, NL-OMON47601). Samples were obtained before treatment (week 0) and at response assessment (week 8). Responders (R) achieved histo-endoscopic mucosal improvement (HEMI) at week 8, defined as an endoscopic Mayo subscore ≤ 1 together with histological remission (Robarts Histopathology Index ≤ 3 without mucosal neutrophils). Patients not achieving HEMI at week 8 were non-responders (NR); the subset of these who continued tofacitinib and showed endoscopic improvement at later routine follow-up were classified as late responders (LR). The transcriptome (bulk RNA-seq of biopsies and whole blood) and the proteome (four Olink Explore 384 panels: Cardiometabolic, Inflammation, Neurology and Oncology; biopsy lysates and serum) were profiled, and the findings were placed in context using two external, publicly available cohorts:

| Modality | Treatment | Source | Platform |
|---|---|---|---|
| Bulk RNA-seq | Tofacitinib | This study | RNA-seq, STAR + featureCounts + DESeq2 |
| Proteomics | Tofacitinib | This study | Olink Explore 384, limma |
| Microarray | Infliximab | Toedter et al. 2011 (GSE23597) | Affymetrix Human Genome U133 Plus 2.0 Array, limma |
| Single-cell RNA-seq (pseudobulk) | Tofacitinib | Melón-Ardanaz et al. 2025 (GSE253006) | 10x Genomics, Seurat pseudobulk + DESeq2 |

Both external datasets differ from the in-house cohort in the clinical metadata available and in their exact definitions of response; these differences are described in the manuscript and should be borne in mind when interpreting cross-cohort comparisons. Agreement across these cohorts is treated as cross-platform corroboration.

---

## Pipeline description

Each modality follows the same conceptual path: prepare a data object, subset per comparison, run differential analysis, run gene-set enrichment:

1. **Bulk RNA-seq (own, tofacitinib).** Copy FASTQ and the STAR index to scratch → FastQC → STAR paired-end alignment on scratch (temporary; copied to the working directory and then removed from scratch) → filter (remove unmapped / multi-mapped / low-MAPQ / mitochondrial reads) with SAMtools → featureCounts (fragment-level) → MultiQC → build a `DESeqDataSet` → per-comparison subsetting → DESeq2 differential expression → fgsea.
2. **Olink proteomics (own, tofacitinib).** Build a `SummarizedExperiment` from NPX values → per-comparison subsetting → limma differential abundance → fgsea.
3. **Microarray (Toedter 2011, infliximab).** Import GSE23597 via GEOquery into a `SummarizedExperiment` → subset → quantile normalization → limma → fgsea.
4. **scRNA-seq pseudobulk (Melón-Ardanaz 2025, tofacitinib).** Pseudobulk each sample with Seurat → merge → build a `DESeqDataSet` → subset → DESeq2 → fgsea.
5. **Interaction test (tofacitinib vs infliximab).** Per-gene Wald test of the difference between the tofacitinib (DESeq2) and infliximab (limma) longitudinal log2 fold changes in responders, with Benjamini-Hochberg correction genome-wide and within the _a priori_ non-canonical NF-κB gene set.
6. **Machine learning.** Pretreatment responder versus non-responder classification per modality (biopsy and blood transcriptome, biopsy and serum proteome) with repeated stratified cross-validation, univariate feature selection and TPOT model tuning confined to the training folds.
7. **Figures and supplementary tables.** Figures 1–5 and Supplementary Tables 1–17 (with legends) are assembled from the outputs above.

**Paired designs.** Longitudinal comparisons are paired on patient. For DESeq2 this is done by including patient as a fixed blocking factor (`~ patient + Week`); patient-level covariates (age, sex, concomitant medication, prior biologics) are constant within patient, absorbed by the blocking factor and therefore dropped automatically, and patients lacking either timepoint are excluded. For limma (Olink, microarray) patient is modelled as a blocking factor via `duplicateCorrelation`. Cross-sectional (baseline) comparisons are adjusted for the listed covariates.

MSigDB gene sets (Hallmark, Reactome and KEGG legacy, via `msigdbr`) are prepared once and shared by every fgsea step.

---

## File structure

Clone the repository and run Snakemake **from the repository root**:

```
<repository root>/
├── workflow/
│   ├── Snakefile                 # entry point; targets defined in `rule all`
│   ├── rules/                    # one .smk module per analysis block
│   │   ├── general.smk           #   MSigDB gene-set preparation (shared)
│   │   ├── rnaseq.smk            #   bulk RNA-seq (tofacitinib)
│   │   ├── olink.smk             #   Olink proteomics (tofacitinib)
│   │   ├── microarray_toedter2011.smk    # infliximab microarray
│   │   ├── scrnaseq_melonardanaz2025.smk # scRNA-seq pseudobulk
│   │   ├── combined.smk          #   tofacitinib vs infliximab interaction test
│   │   ├── machine_learning.smk  #   pretreatment response prediction
│   │   ├── figures.smk           #   figures 1–5
│   │   └── supplementary_tables.smk      # supplementary tables 1–17 + legends
│   ├── scripts/                  # R / Python analysis scripts, one folder per block
│   └── envs/                     # per-rule conda environment definitions
│
├── config/
│   ├── config.yaml               # all paths, parameters and comparisons
│   ├── supplementary_tables.yaml # supplementary table definitions and legends
│   ├── metadata/                 # file manifests (xlsx); sample metadata not included
│   └── features_of_interest/     # curated gene and gene-set lists (xlsx)
│
├── resources/                    # external raw data — you provide these
│   ├── toedter2011/              #   GSE23597 series matrix
│   └── melonardanaz2025/         #   GSE253006 matrices, one dir per sample
│
├── output/                       # created by the workflow
└── supplementary_tables/         # created by the workflow
```

---

## Software requirements

- **Snakemake** (workflow manager) with the **conda** integration enabled (`--software-deployment-method conda`). Each rule pins its own environment, so you do not need to install the scientific stack yourself. Snakemake creates the environments on first run.
- **conda / mamba** to solve those environments.
- A **SLURM** cluster is assumed for the sequencing-heavy steps (the study was run on a BeeGFS-backed HPC running RHEL). The workflow runs locally too, but STAR alignment, pseudobulk merging and the machine-learning sweep have substantial memory or compute requirements (see the `resources:` directives in the rules).

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
| Python | 3.10 |
| scikit-learn | 1.3.2 |
| TPOT | 0.12.2 |
| XGBoost | 3.0.5 |
| pandas / openpyxl (supplementary tables) | 2.2.3 / 3.1.5 |

---

## Input data

The inputs fall into three categories, and this matters for reproducibility:

**1. Provided in this repository** — file manifests (`config/metadata/rnaseq/`, `config/metadata/olink/`, `config/metadata/scrnaseq_melonardanaz2025/`) and curated gene/gene-set lists (`config/features_of_interest/`).

**Sample metadata (not included)** — the per-sample clinical metadata used for the comparisons and covariates (`data.sample_metadata_xlsx` in `config.yaml`, expected at `config/metadata/samples/`) contains pseudonymised patient-level information and is therefore not distributed in this repository. It is available on reasonable request (see [License and contact](#license-and-contact)).

**2. Public, must be downloaded into `resources/`:**

- **GSE23597** (Toedter et al. 2011, infliximab microarray). Place the series matrix at the path given by `data.toedter2011_gse` in `config.yaml`
  (`resources/toedter2011/GSE23597_series_matrix.txt`).
- **GSE253006** (Melón-Ardanaz et al. 2025, scRNA-seq). Place the per-sample count matrices under `data.melonardanaz2025_mtx_dir`
  (`resources/melonardanaz2025/<Sample_ID>/`).

**3. Not included** — the in-house raw sequencing reads (FASTQ) and Olink NPX measurements are individual-level human data and are not distributed here. Processed RNA-seq counts, Olink NPX values and all summary statistics are deposited on Zenodo (DOI: [10.5281/zenodo.23210920](https://doi.org/10.5281/zenodo.23210920)). Without the raw data, the RNA-seq alignment and Olink preparation steps cannot be re-run, but every downstream differential-expression and enrichment step can be reproduced from the processed data.

**Reference genome/annotation.** The STAR index (GRCh38, Ensembl 104) and GTF (GENCODE v38) are referenced in `config.yaml` by absolute cluster paths under `reference:`. Point these at your own local copies before running the alignment rules.

---

## Configuration

Everything the workflow needs is declared in `config/config.yaml`:

- `base_dir`, `scratch_dir` — project root and fast scratch location on the HPC. Edit these for your environment.
- `R:` — shared plotting parameters and the curated gene / gene-set spreadsheets used by the figure scripts.
- `data:` — paths to every metadata manifest and external dataset.
- `reference:` — STAR index and GTF paths.
- `analysis_params:` — filtering and enrichment thresholds (DESeq2 count filter, limma detection filter, fgsea min/max set size and permutations).
- `interaction:` — the _a priori_ gene set used by the tofacitinib-vs-infliximab interaction test, and the comparisons it is run for.
- `machine_learning:` — modalities, models, feature-set sizes, number of repeats and TPOT settings for the response-prediction sweep.
- Three comparison blocks (see below).

Supplementary tables (source output, comparison and legend per table) are defined in `config/supplementary_tables.yaml`.

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
| `block` | _(optional)_ within-subject blocking factor (`DonorID`, `Donor_ID` for GSE253006) for paired designs. Only used when the tested variable varies within the block. |
| `covariates` | _(optional)_ `numeric:` and `categorical:` adjustment covariates. In paired DESeq2 designs, covariates that are constant within the block are dropped. |

---

## Comparison naming convention

Comparison identifiers encode the compartment, the response group and the contrast:

```
  BINF_R_W8vW0
  │    │ └── contrast: week 8 vs week 0 (longitudinal)
  │    └──── response group
  └───────── compartment
```

**Compartment:** `BINF` = inflamed biopsy. `BNINF` = non-inflamed biopsy. `PB` = peripheral blood (whole-blood RNA-seq; serum for Olink).

**Response group:** `R` = responder. `NR` = non-responder. `LR` = late responder.

**Contrast:** `W8vW0` = longitudinal (week 8 vs week 0 within one response group). `RvNR_W0` = cross-sectional (responders vs non-responders at baseline).

---

## Running the workflow

From the repository root, with conda/mamba and Snakemake available:

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

The default targets are declared in `rule all` in the `Snakefile`. Many intermediate targets are commented out there so that a full run produces the enrichment tables, the interaction results, the machine-learning summaries, the figures and the supplementary tables. Uncomment the relevant lines (or name a specific output on the command line) to materialize intermediate objects such as BAMs, count matrices or per-comparison `DESeqDataSet`/`SummarizedExperiment` objects.

---

## Outputs

All results are written under `output/`, organized by modality and comparison, for example:

```
output/
├── general/gs_msigdb.Rds
├── rnaseq/
│   ├── bam/<sample>/            (STAR alignments, copied from scratch)
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
├── combined/                    (tofacitinib vs infliximab interaction tables)
├── machine_learning/{data,raw,feature_selection,summary}/
└── figures/figure{1..5}/*.pdf

supplementary_tables/
├── tableS1.xlsx … tableS17.xlsx
└── supplementary_table_legends.txt
```

Differential tables carry effect sizes, standard errors and adjusted p-values. Fgsea tables carry normalized enrichment scores, adjusted p-values and pipe-delimited leading-edge gene lists.

---

## Notes, limitations and caveats

- **External-cohort differences.** The infliximab (Toedter 2011) and scRNA-seq (Melón-Ardanaz 2025) datasets lack some clinical covariates present in the in-house cohort and define response differently. Adjustment models for these cohorts are therefore simpler than for the in-house data, and cross-cohort comparisons are interpreted with that asymmetry in mind.
- **Cross-platform interaction test.** Tofacitinib effects are estimated on RNA-seq and infliximab effects on microarray, whose log2 fold-change scales are not identical; the interaction script also reports a quantile-normalized effect difference as a scale-robustness check.
- **Machine learning.** Performance is estimated by repeated cross-validation within a small cohort without an external validation set, and the reported configuration per modality is the best of several models and feature-set sizes. Reported AUROCs are therefore optimistic and should be treated as a relative comparison of modalities rather than as expected performance in new patients.
- **Reproducibility scope.** Steps upstream of the differential-expression tables depend on individual-level raw data (see [Input data](#input-data)). The enrichment, interaction, figure and supplementary-table steps are reproducible from the processed data.

---

## Generative AI usage statement

Parts of this project were developed with the assistance of Anthropic Claude, Opus 4.8. It was used for generating boilerplate code, debugging, and drafting documentation (including this README.md). All code and documentation has been reviewed, further developed, and tested by a human.

---

## Citation

If you use this code, please cite the manuscript:

```
Fung ICN, de Jong DC, Liu K, et al. Longitudinal multi-omics profiling of tofacitinib response identifies a pre-existing Janus kinase-tractable inflammatory state in ulcerative colitis. medRxiv 2026. DOI: 10.64898/2026.10.07.26364964
[TODO: replace with the journal citation on acceptance]
```

and the external datasets it builds on:

- Toedter G, et al. 2011. Gene expression profiling and response signatures associated with differential responses to infliximab treatment in ulcerative colitis. _Am J Gastroenterol_ 2011. DOI: 10.1038/ajg.2011.83. GEO Accession ID: GSE23597.
- Melón-Ardanaz E, et al. 2025. Differential effects of tofacitinib on macrophage activation contribute to lack of response in ulcerative colitis patients. _J Crohns Colitis_. DOI: 10.1093/ecco-jcc/jjaf076. GEO Accession ID: GSE253006.

---

## License and contact

- **License:** MIT (see `LICENSE`)
- **Corresponding author / maintainer:** Andrew Y.F. Li Yim, ORCID: 0000-0002-0754-0953, Email: a.y.liyim [at] amsterdamumc.nl
