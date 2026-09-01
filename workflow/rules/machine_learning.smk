# =============================================================================
# machine_learning.smk  —  pretreatment response prediction (modality comparison)
# -----------------------------------------------------------------------------
# The predictive sweep is expressed natively as a Snakemake job array: one
# `ml_train` instance per (modality x model x k). Run under the SLURM executor,
# Snakemake submits and tracks these directly (and `--rerun-incomplete`
# resubmits failures), which replaces the standalone submit_jobs.sh /
# monitor_and_resubmit.sh orchestration.
#
#   ml_prepare_data   dds.Rds + se.Rds        -> feature matrices + labels.csv
#   ml_train          feature matrix + labels -> per-(modality,model,k) folds
#   ml_aggregate      all fold CSVs           -> summary tables + plots
# =============================================================================

# modality -> exported feature-matrix filename (written by ml_prepare_data)
ML_DATA_FILE = {
    "Tissue_RNA": "rna_tissue_data.csv",
    "PBL_RNA": "rna_pbl_data.csv",
    "Tissue_PROT": "prot_tissue_data.csv",
    "PBL_PROT":    "prot_pbl_data.csv",
}

rule ml_prepare_data:
  input:
    dds_rds="output/rnaseq/dds/dds.Rds",
    se_rds="output/olink/se/se.Rds",
  output:
    labels_csv="output/machine_learning/data/RvNR_W0_labels.csv",
    tx_binf_rvnr_w0_data_csv="output/machine_learning/data/TX_BINF_RvNR_W0_data.csv",
    tx_pb_rvnr_w0_data_csv="output/machine_learning/data/TX_PB_RvNR_W0_data.csv",
    px_binf_rvnr_w0_data_csv="output/machine_learning/data/PX_BINF_RvNR_W0_data.csv",
    px_pb_rvnr_w0_data_csv="output/machine_learning/data/PX_PB_RvNR_W0_data.csv",
  params:
    ml_params=config["machine_learning"],
  log:
    "output/machine_learning/data/ml_prepare_data.log",
  conda:
    "../envs/r-deseq2.yaml"
  message:
    "--- Preparing ML feature data matrices and labels ---"
  threads: 1
  resources:
    mem_mb=16000,
  script:
    "../scripts/machine_learning/prepare_data.R"

rule ml_train:
  input:
    data=lambda wc: f"output/machine_learning/data/{ML_DATA_FILE[wc.modality]}",
    labels="output/machine_learning/data/labels.csv",
  output:
    folds_csv="output/machine_learning/raw/{modality}__{model}__k{k}.csv",
  params:
    modality=lambda wc: wc.modality,
    model=lambda wc: wc.model,
    k=lambda wc: int(wc.k),
    metric=_ml["metric"],
    n_seeds=_ml["n_seeds"],
    seed=_ml["seed"],
    tpot=_ml["tpot"],
  log:
    "output/machine_learning/raw/logs/{modality}__{model}__k{k}.log",
  conda:
    "../envs/machine_learning.yaml"
  message:
    "--- ML train: {wildcards.modality} / {wildcards.model} / k={wildcards.k} ---"
  threads: 2
  resources:
    mem_mb=16000,
    runtime=1440,   # minutes; each job runs n_seeds x LOOCV x TPOT
  script:
    "../scripts/machine_learning/worker.py"

rule ml_aggregate:
  input:
    folds=expand(
      "output/machine_learning/raw/{modality}__{model}__k{k}.csv",
      modality=_ml["modalities"], model=_ml["models"], k=_ml["k_values"],
    ),
  output:
    fold_level="output/machine_learning/aggregated_fold_level.csv",
    summary_by_k="output/machine_learning/summary_by_k.csv",
    summary_best_k="output/machine_learning/summary_best_k.csv",
    roc_plot="output/machine_learning/plots/roc_curves_combined.png",
  params:
    ml_params=_ml,
  log:
    "output/machine_learning/ml_aggregate.log",
  conda:
    "../envs/machine_learning.yaml"
  message:
    "--- Aggregating ML results ---"
  threads: 1
  resources:
    mem_mb=16000,
  script:
    "../scripts/machine_learning/aggregate_results.py"
