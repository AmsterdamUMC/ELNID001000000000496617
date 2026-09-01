#!/usr/bin/env bash
# =============================================================================
# monitor_and_resubmit.sh
# Check job completion and resubmit any failed/missing jobs.
#
# Usage:
#   bash monitor_and_resubmit.sh          # show status only
#   bash monitor_and_resubmit.sh --fix    # resubmit missing/failed jobs
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
RESULT_DIR="${BASE_DIR}/results/raw"
SUBMIT_SCRIPT="${SCRIPT_DIR}/submit_jobs.sh"

RESUBMIT=false
for arg in "$@"; do
    [[ "$arg" == "--fix" ]] && RESUBMIT=true
done

MODALITIES=("PBL_PROT" "Tissue_PROT" "PBL_RNA" "Tissue_RNA")
MODEL_NAMES=("LogisticRegression" "GaussianNB" "XGBoost" "KNeighbors")
K_VALUES=(10 25 50 100 200 300 400 500)

EXPECTED=$(( ${#MODALITIES[@]} * ${#MODEL_NAMES[@]} * ${#K_VALUES[@]} ))
DONE=0
MISSING=0

for MOD in "${MODALITIES[@]}"; do
    for MNAME in "${MODEL_NAMES[@]}"; do
        for K in "${K_VALUES[@]}"; do
            OUT_CSV="${RESULT_DIR}/${MOD}__${MNAME}__k${K}.csv"
            if [[ -f "$OUT_CSV" ]]; then
                (( DONE++ )) || true
            else
                (( MISSING++ )) || true
            fi
        done
    done
done

echo "══════════════════════════════════════"
echo " Expected  : ${EXPECTED}"
echo " Done (OK) : ${DONE}"
echo " Missing   : ${MISSING}"
echo " Progress  : $(echo "scale=1; ${DONE} * 100 / ${EXPECTED}" | bc)%"
echo "══════════════════════════════════════"

echo ""
echo "SLURM queue (your jobs):"
squeue --me --format="%.10i %.20j %.8T %.10M %.5C" 2>/dev/null || echo "(squeue not available)"

if [[ "$RESUBMIT" == true && $MISSING -gt 0 ]]; then
    echo ""
    echo "Resubmitting ${MISSING} missing jobs..."
    bash "$SUBMIT_SCRIPT" --resume
fi

if [[ $MISSING -eq 0 ]]; then
    echo ""
    echo "✅ All jobs complete! Run:"
    echo "   python ${SCRIPT_DIR}/aggregate_results.py"
fi
