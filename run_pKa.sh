#!/usr/bin/env bash
set -euo pipefail

# Get script directory (safe even if run from elsewhere)
SCRIPT_DIR="$(cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P)"
RESULTS_DIR="$SCRIPT_DIR/results/lead_generation_output"

# Input/Output
IN="$RESULTS_DIR/generated_smiles_pass_pH6_3.smi"    # input: result filtered at pH 6.3
OUT="$RESULTS_DIR/generated_smiles_pass_pKa.smi"     # final output after pKa filter
LOG="$RESULTS_DIR/pka_log.csv"

# Create results directory if missing
mkdir -p "$RESULTS_DIR"

# ---- pKa filtering criteria ----
PKA_MIN=5
PKA_MAX=8
CKPT="./datasets/data/pKa/last.ckpt"
# --------------------------------

# Run pKa filter
python utils/predict_pka.py \
  --in-smi "$IN" \
  --out-csv "$LOG" \
  --out-smi "$OUT" \
  --reg-ckpt "$CKPT" \
  --reg-min "$PKA_MIN" \
  --reg-max "$PKA_MAX"

echo "[pKa filter] Final passed: $(wc -l < "$OUT" | tr -d ' ') molecules -> $OUT"
