comparisons_olink = [c for c in comparisons if 'olink'  in comparisons[c].get('assays', ['rnaseq', 'olink'])]

# Rules

rule olink_se_preparation:
  input:
    olink_files_xlsx=config["data"]["olink_files_xlsx"],
    sample_metadata_xlsx=config["data"]["sample_metadata_xlsx"],
  output:
    se_rds="output/olink/se/se.Rds",
  conda:
    "../envs/r-summarizedexperiment.yaml"
  log:
    "output/olink/se/olink_se_preparation.log",
  message:
    "--- Preparing SummarizedExperiment OLINK ---"
  threads: 1
  resources:
    mem_mb=8000,
  script:
    "../scripts/olink/prepare_se.R"

rule olink_se_subset:
  input:
    se_rds="output/olink/se/se.Rds",
  output:
    se_subset_rds="output/olink/analyses/{comparison}/se/se_{comparison}.Rds",
  params:
    filter=lambda wc: config["comparisons"][wc.comparison]["filter"]
  log:
    "output/olink/analyses/{comparison}/se/olink_se_{comparison}_subsetting.log",
  conda:
    "../envs/r-summarizedexperiment.yaml"
  message:
    "--- Subsetting {wildcards.comparison} OLINK ---"
  threads: 1
  resources:
    mem_mb=16000,
  script:
    "../scripts/olink/subset_se.R"

rule olink_differential_expression:
  input:
    se_rds="output/olink/analyses/{comparison}/se/se_{comparison}.Rds",
  output:
    marraylm_rds="output/olink/analyses/{comparison}/deps/marraylm_{comparison}.Rds",
    deps_csv="output/olink/analyses/{comparison}/deps/deps_{comparison}.csv",
  params:
    comparison=lambda wc: config["comparisons"][wc.comparison],
    limma_params=lambda wc: config["analysis_params"]["limma"],
  log:
    "output/olink/analyses/{comparison}/deps/olink_differential_expression_{comparison}.log",
  conda:
    "../envs/r-limma.yaml"
  message:
    "--- Differential expression {wildcards.comparison} OLINK ---"
  threads: 1
  resources:
    mem_mb=16000,
  script:
    "../scripts/olink/differential_expression.R"

rule olink_gene_set_enrichment_analysis:
  input:
    msigdb_rds="output/general/gs_msigdb.Rds",
    deps_csv="output/olink/analyses/{comparison}/deps/deps_{comparison}.csv",
  output:
    fgsea_csv="output/olink/analyses/{comparison}/fgsea/fgsea_{comparison}.csv",
  params:
    fgsea_params=lambda wc: config["analysis_params"]["fgsea"],
  log:
    "output/olink/analyses/{comparison}/fgsea/olink_gene_set_enrichment_analysis_{comparison}.log",
  conda:
    "../envs/r-fgsea.yaml"
  message:
    "--- Gene set enrichment analysis {wildcards.comparison} OLINK ---"
  threads: 1
  resources:
    mem_mb=16000,
  script:
    "../scripts/olink/fgsea.R"
