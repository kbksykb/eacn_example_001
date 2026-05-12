#!/bin/bash
# launch_kang2024_fig4.sh — launch all NMF compartments in parallel, one per GPU.
# Usage: bash launch_kang2024_fig4.sh
# Logs go to workspace/logs/kang2024_fig4_<compartment>.log

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="$SCRIPT_DIR/kang2024_fig4_pipeline.py"
LOG_DIR="$SCRIPT_DIR/../../logs"
mkdir -p "$LOG_DIR"

CONDA_PYTHON="/root/anaconda3/envs/scenv/bin/python"

declare -A COMPARTMENT_GPU=(
  ["tnk"]="0"
  ["myl"]="1"
  ["b"]="2"
  ["mesenchymal"]="3"
)

PIDS=()
for compartment in "${!COMPARTMENT_GPU[@]}"; do
  gpu="${COMPARTMENT_GPU[$compartment]}"
  logfile="$LOG_DIR/kang2024_fig4_${compartment}.log"
  echo "Launching $compartment on GPU $gpu -> $logfile"
  CUDA_VISIBLE_DEVICES=$gpu $CONDA_PYTHON "$PIPELINE" \
    --compartment "$compartment" \
    > "$logfile" 2>&1 &
  PIDS+=($!)
done

echo "All compartments launched. PIDs: ${PIDS[*]}"
echo "Monitor with: tail -f $LOG_DIR/kang2024_fig4_*.log"

# Wait for all
for pid in "${PIDS[@]}"; do
  wait "$pid" && echo "PID $pid done OK" || echo "PID $pid FAILED"
done
echo "All compartments complete."
