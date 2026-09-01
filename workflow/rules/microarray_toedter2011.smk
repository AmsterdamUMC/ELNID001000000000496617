# Toedter et al. 2011 (GSE23597) pertains a bulk transcriptomic microarray study that performed a longitudinal case-control study on inflamed intestinal biopsies obtained from ulcerative colitis patients starting anti-TNF treatment infliximab treatment. Samples were obtained at week 0 (pretreatment), 8, and 30 from patients treated with either infliximab or placebo. We note that the GSE does not encompass certain clinical parameters that we included in our own study, such as age, sex, and/or prior anti-TNF or vedolizumab usage. I doubt that any were on vedolizumab as that was not approved back then. It is possible some were adalimumab-exposed, but no information is available. Another point of concern is their definitions of response, which differ slightly from our own, given that they consider total Mayo score as a compound of its constituents: we cannot disentangle the individual consituents.

comparisons_toedter2011 = config['comparisons_toedter2011']

# Rules

rule microarray_toedter2011_se_preparation:
  input:
    toedter2011_gse=config["data"]["toedter2011_gse"],
  output:
    se_rds="output/microarray_toedter2011/se/se.Rds",
  conda:
    "../envs/r-geoquery.yaml"
  log:
    "output/microarray_toedter2011/se/microarray_toedter2011_infliximab_se_preparation.log",
  message:
    "--- Preparing SummarizedExperiment Toedter 2011 ---"
  threads: 1
  resources:
    mem_mb=8000,
  script:
    "../scripts/microarray_toedter2011/prepare_se.R"

rule microarray_toedter2011_se_subset:
  input:
    se_rds="output/microarray_toedter2011/se/se.Rds",
  output:
    se_subset_rds="output/microarray_toedter2011/analyses/{comparison}/se/se_{comparison}.Rds",
  params:
    filter=lambda wc: config["comparisons_toedter2011"][wc.comparison]["filter"]
  log:
    "output/microarray_toedter2011/analyses/{comparison}/se/microarray_toedter2011_infliximab_se_{comparison}_subsetting.log",
  conda:
    "../envs/r-limma.yaml"
  message:
    "--- Subsetting {wildcards.comparison} Toedter 2011 ---"
  threads: 1
  resources:
    mem_mb=16000,
  script:
    "../scripts/microarray_toedter2011/subset_se.R"

rule microarray_toedter2011_differential_expression:
  input:
    se_rds="output/microarray_toedter2011/analyses/{comparison}/se/se_{comparison}.Rds",
  output:
    marraylm_rds="output/microarray_toedter2011/analyses/{comparison}/degs/marraylm_{comparison}.Rds",
    degs_csv="output/microarray_toedter2011/analyses/{comparison}/degs/degs_{comparison}.csv",
  params:
    comparison=lambda wc: config["comparisons_toedter2011"][wc.comparison],
  log:
    "output/microarray_toedter2011/analyses/{comparison}/degs/microarray_toedter2011_differential_expression_{comparison}.log",
  conda:
    "../envs/r-limma.yaml"
  message:
    "--- Differential expression {wildcards.comparison} Toedter 2011 ---"
  threads: 1
  resources:
    mem_mb=16000,
  script:
    "../scripts/microarray_toedter2011/differential_expression.R"

rule microarray_toedter2011_gene_set_enrichment_analysis:
  input:
    msigdb_rds="output/general/gs_msigdb.Rds",
    degs_csv="output/microarray_toedter2011/analyses/{comparison}/degs/degs_{comparison}.csv",
  output:
    fgsea_csv="output/microarray_toedter2011/analyses/{comparison}/fgsea/fgsea_{comparison}.csv",
  params:
    fgsea_params=lambda wc: config["analysis_params"]["fgsea"],
  log:
    "output/microarray_toedter2011/analyses/{comparison}/fgsea/microarray_toedter2011_gene_set_enrichment_analysis_{comparison}.log",
  conda:
    "../envs/r-fgsea.yaml"
  message:
    "--- Gene set enrichment analysis {wildcards.comparison} Toedter 2011 ---"
  threads: 1
  resources:
    mem_mb=16000,
  script:
    "../scripts/microarray_toedter2011/fgsea.R"
