#!/usr/bin/env bash
# =============================================================================
# submit_jobs.sh
# Submit one SLURM job per (modality × model × k).
#
# Total jobs = 4 modalities × 4 models × 8 k-values = 128 jobs
# Each job runs 100 seeds × LOOCV internally (~8–24 h per job).
#
# Usage:
#   bash submit_jobs.sh              # submit all
#   bash submit_jobs.sh --dry-run    # print commands, don't submit
#   bash submit_jobs.sh --resume     # skip already-completed results
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
WORKER="${SCRIPT_DIR}/worker.py"
LOG_DIR="${BASE_DIR}/logs"
RESULT_DIR="${BASE_DIR}/results/raw"

mkdir -p "$LOG_DIR" "$RESULT_DIR"

# ── SLURM resource settings — tune to your cluster ────────────
PARTITION="defq"           # 100 seeds × LOOCV needs a long partition
TIME="24:00:00"            # wall time per job
MEM="16G"
CPUS=2

# ── Sweep parameters ──────────────────────────────────────────
MODALITIES=("PBL_PROT" "Tissue_PROT" "PBL_RNA" "Tissue_RNA")
MODEL_IDXS=(0 1 2 3)
MODEL_NAMES=("LogisticRegression" "GaussianNB" "XGBoost" "KNeighbors")
K_VALUES=(10 25 50 100 200 300 400 500)

# ── Flags ─────────────────────────────────────────────────────
DRY_RUN=false
RESUME=false
for arg in "$@"; do
    [[ "$arg" == "--dry-run" ]] && DRY_RUN=true
    [[ "$arg" == "--resume"  ]] && RESUME=true
done

# ── Submit ────────────────────────────────────────────────────
n_submitted=0
n_skipped=0

for MOD in "${MODALITIES[@]}"; do
    for i in "${!MODEL_IDXS[@]}"; do
        MIDX="${MODEL_IDXS[$i]}"
        MNAME="${MODEL_NAMES[$i]}"
        for K in "${K_VALUES[@]}"; do

            OUT_CSV="${RESULT_DIR}/${MOD}__${MNAME}__k${K}.csv"

            if [[ "$RESUME" == true && -f "$OUT_CSV" ]]; then
                (( n_skipped++ )) || true
                continue
            fi

            JOB_NAME="${MOD}_m${MIDX}_k${K}"
            LOG_FILE="${LOG_DIR}/${JOB_NAME}_%j.log"

            SBATCH_CMD=(
                sbatch
                --job-name="$JOB_NAME"
                --partition="$PARTITION"
                --time="$TIME"
                --mem="$MEM"
                --cpus-per-task="$CPUS"
                --output="$LOG_FILE"
                --wrap=" source /appdata/users/P004033/miniforge3/etc/profile.d/mamba.sh
                    mamba activate AI && \
                    python ${WORKER} \
                        --modality  ${MOD}  \
                        --model_idx ${MIDX} \
                        --k         ${K}    \
                        --metric    roc_auc"
            )

            if [[ "$DRY_RUN" == true ]]; then
                echo "[DRY-RUN] ${SBATCH_CMD[*]}"
            else
                "${SBATCH_CMD[@]}"
            fi

            (( n_submitted++ )) || true
        done
    done
done

echo ""
echo "──────────────────────────────────────────"
echo " Jobs submitted : ${n_submitted}"
echo " Jobs skipped   : ${n_skipped}  (already done)"
echo " Logs           : ${LOG_DIR}"
echo " Results        : ${RESULT_DIR}"
echo "──────────────────────────────────────────"
