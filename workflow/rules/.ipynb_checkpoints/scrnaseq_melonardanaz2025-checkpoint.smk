# Melon-Ardanaz et al. 2011 (GSE253006) pertains a single-cell RNA-sequencing study that performed a longitudinal case-control study on inflamed intestinal biopsies obtained from ulcerative colitis patients starting tofacitinib treatment. Samples were obtained at week 0 (pretreatment) and 8 (with some at week 24 and 48) from patients (starting treatment) with tofacitinib. We note that the GSE does not encompass certain clinical parameters that we included in our own study, such as age, sex, and/or prior anti-TNF or vedolizumab usage. Another point of concern is their definitions of response, which differ slightly from our own, where they define response as a decrease in eMayo =< 1 as response.

comparisons_scrnaseq_melonardanaz2025 = config['comparisons_melonardanaz2025']
scrnaseq_melonardanaz2025_files_df = pd.read_excel(config['data']['melonardanaz2025_files_xlsx'], header=0)
scrnaseq_melonardanaz2025_sampleIDs = pd.Series(scrnaseq_melonardanaz2025_files_df['Sample_ID'].unique())

rule scrnaseq_melonardanaz2025_seuratobject_preparation:
  input:
    melonardanaz2025_dir=config["data"]["melonardanaz2025_mtx_dir"] + "/{scrnaseq_melonardanaz2025_sampleID}",
  output:
    pbcounts_csv="output/scrnaseq_melonardanaz2025/import/{scrnaseq_melonardanaz2025_sampleID}_pseudobulk.csv",
  threads:
    1
  conda:
    "../envs/r-seurat.yaml"
  log:
    "output/scrnaseq_melonardanaz2025/import/scrnaseq_melonardanaz2025_{scrnaseq_melonardanaz2025_sampleID}_preparation.log",
  benchmark:
    "output/scrnaseq_melonardanaz2025/import/scrnaseq_melonardanaz2025_{scrnaseq_melonardanaz2025_sampleID}_preparation_benchmark.txt"
  resources:
    mem_mb=16000,
  params:
    sampleid="{scrnaseq_melonardanaz2025_sampleID}",
    scratch_dir=config["scratch_dir"],
  script:
    "../scripts/scrnaseq_melonardanaz2025/prepare_pseudobulk.R"

rule scrnaseq_melonardanaz2025_seuratobject_merge:
  input:
    pbcounts_csvs=expand("output/scrnaseq_melonardanaz2025/import/{scrnaseq_melonardanaz2025_sampleID}_pseudobulk.csv", scrnaseq_melonardanaz2025_sampleID=scrnaseq_melonardanaz2025_sampleIDs),
  output:
    pbcounts_merged_csv="output/scrnaseq_melonardanaz2025/merged/pseudobulk_merged.csv",
  threads:
    1
  conda:
    "../envs/r-seurat.yaml",
  log:
    "output/scrnaseq_melonardanaz2025/merged/scrnaseq_melonardanaz2025_seuratobject_merge.log",
  benchmark:
    "output/scrnaseq_melonardanaz2025/merged/scrnaseq_melonardanaz2025_seuratobject_merge_benchmark.txt",
  resources:
    mem_mb=300000,
  script:
    "../scripts/scrnaseq_melonardanaz2025/merge_pseudobulk.R"

rule scrnaseq_melonardanaz2025_dds_preparation:
  input:
    counts_csv="output/scrnaseq_melonardanaz2025/merged/pseudobulk_merged.csv",
    sample_metadata_xlsx=config["data"]["melonardanaz2025_metadata_xlsx"],
  output:
    dds_rds="output/scrnaseq_melonardanaz2025/dds/dds.Rds",
  conda:
    "../envs/r-deseq2.yaml"
  log:
    "output/scrnaseq_melonardanaz2025/dds/rnaseq_dds_preparation.log",
  message:
    "--- Preparing DESeqDataSet scRNAseq ---"
  threads: 1
  resources:
    mem_mb=8000,
  script:
    "../scripts/scrnaseq_melonardanaz2025/prepare_dds.R"

rule scrnaseq_melonardanaz2025_dds_subset:
  input:
    dds_rds="output/scrnaseq_melonardanaz2025/dds/dds.Rds",
  output:
    dds_subset_rds="output/scrnaseq_melonardanaz2025/analyses/{comparison_scrnaseq_melonardanaz2025}/dds/dds_{comparison_scrnaseq_melonardanaz2025}.Rds",
    vst_subset_rds="output/scrnaseq_melonardanaz2025/analyses/{comparison_scrnaseq_melonardanaz2025}/dds/vst_{comparison_scrnaseq_melonardanaz2025}.Rds",
  params:
    filter=lambda wc: config["comparisons_melonardanaz2025"][wc.comparison_scrnaseq_melonardanaz2025]["filter"]
  log:
    "output/scrnaseq_melonardanaz2025/analyses/{comparison_scrnaseq_melonardanaz2025}/dds/rnaseq_dds_{comparison_scrnaseq_melonardanaz2025}_subsetting.log",
  conda:
    "../envs/r-deseq2.yaml"
  message:
    "--- Subsetting {wildcards.comparison_scrnaseq_melonardanaz2025} scRNAseq ---"
  threads: 1
  resources:
    mem_mb=16000,
  script:
    "../scripts/scrnaseq_melonardanaz2025/subset_dds.R"

rule scrnaseq_melonardanaz2025_differential_expression:
  input:
    dds_rds="output/scrnaseq_melonardanaz2025/analyses/{comparison_scrnaseq_melonardanaz2025}/dds/dds_{comparison_scrnaseq_melonardanaz2025}.Rds",
  output:
    deseq_rds="output/scrnaseq_melonardanaz2025/analyses/{comparison_scrnaseq_melonardanaz2025}/degs/deseq_{comparison_scrnaseq_melonardanaz2025}.Rds",
    degs_csv="output/scrnaseq_melonardanaz2025/analyses/{comparison_scrnaseq_melonardanaz2025}/degs/degs_{comparison_scrnaseq_melonardanaz2025}.csv",
  params:
    comparison=lambda wc: config["comparisons_melonardanaz2025"][wc.comparison_scrnaseq_melonardanaz2025],
    deseq2_params=lambda wc: config["analysis_params"]["deseq2"],
  log:
    "output/scrnaseq_melonardanaz2025/analyses/{comparison_scrnaseq_melonardanaz2025}/degs/rnaseq_differential_expression_{comparison_scrnaseq_melonardanaz2025}.log",
  conda:
    "../envs/r-deseq2.yaml"
  message:
    "--- Differential expression {wildcards.comparison_scrnaseq_melonardanaz2025} scRNAseq ---"
  threads: 1
  resources:
    mem_mb=16000,
  script:
    "../scripts/scrnaseq_melonardanaz2025/differential_expression.R"

rule scrnaseq_melonardanaz2025_gene_set_enrichment_analysis:
  input:
    msigdb_rds="output/general/gs_msigdb.Rds",
    degs_csv="output/scrnaseq_melonardanaz2025/analyses/{comparison_scrnaseq_melonardanaz2025}/degs/degs_{comparison_scrnaseq_melonardanaz2025}.csv",
  output:
    fgsea_csv="output/scrnaseq_melonardanaz2025/analyses/{comparison_scrnaseq_melonardanaz2025}/fgsea/fgsea_{comparison_scrnaseq_melonardanaz2025}.csv",
  params:
    fgsea_params=lambda wc: config["analysis_params"]["fgsea"],
  log:
    "output/scrnaseq_melonardanaz2025/analyses/{comparison_scrnaseq_melonardanaz2025}/fgsea/rnaseq_gene_set_enrichment_analysis_{comparison_scrnaseq_melonardanaz2025}.log",
  conda:
    "../envs/r-fgsea.yaml"
  message:
    "--- Gene set enrichment analysis {wildcards.comparison_scrnaseq_melonardanaz2025} scRNAseq ---"
  threads: 1
  resources:
    mem_mb=16000,
  script:
    "../scripts/scrnaseq_melonardanaz2025/fgsea.R"
