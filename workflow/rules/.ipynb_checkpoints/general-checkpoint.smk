# Rules

rule prepare_msigdbr:
  output:
    gs_rds="output/general/gs_msigdb.Rds",
  log:
    "output/general/prepare_msigdbr.log",
  conda:
    "../envs/r-msigdbr.yaml"
  message:
    "--- Download MSigDBR gene sets ---"
  threads: 1
  resources:
    mem_mb=16000,
  script:
    "../scripts/prepare_msigdbr.R"
