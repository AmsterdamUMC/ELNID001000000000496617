rule interaction_infliximab_tofacitinib_analysis:
  input:
    infliximab_degs_csv="output/microarray_toedter2011/analyses/{comparison_interaction}/degs/degs_{comparison_interaction}.csv",
    tofacitinib_degs_csv="output/rnaseq/analyses/{comparison_interaction}/degs/degs_{comparison_interaction}.csv",
  output:
    interaction_degs_csv="output/combined/interaction_infliximab_tofacitinib_{comparison_interaction}.csv",
    interaction_nfkbnc_csv="output/combined/interaction_infliximab_tofacitinib_nfkbnc_{comparison_interaction}.csv",
  log:
    "output/combined/infliximab_tofacitinib_interaction_analysis_{comparison_interaction}.log",
  conda:
    "../envs/r-deseq2.yaml"
  message:
    "--- Interaction: tofacitinib vs infliximab over responders over time analysis {wildcards.comparison_interaction} ---"
  threads: 1
  resources:
    mem_mb=16000,
  script:
    "../scripts/combined/interaction_tofacitinib_infliximab.R"
