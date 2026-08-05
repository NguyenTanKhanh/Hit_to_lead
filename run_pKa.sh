#!/usr/bin/env bash
set -euo pipefail

# ===============================
# Hardening: avoid AMBER/MPI hijack
# ===============================
unset PYTHONPATH
unset PYTHONHOME
unset PYTHONUSERBASE
export PYTHONNOUSERSITE=1

# Disable Lightning environment detection (MPI/SLURM)
export PL_DISABLE_ENVIRONMENT_DETECTION=1

# Unset common MPI/SLURM variables that trigger Lightning detection
unset SLURM_JOB_ID
unset SLURM_PROCID
unset PMI_RANK
unset PMIX_RANK
unset OMPI_COMM_WORLD_RANK
unset OMPI_COMM_WORLD_SIZE
unset MPI_LOCALRANKID
unset MPI_LOCALNRANKS
unset I_MPI_ROOT
unset I_MPI_INFO

# Optional: make sure we are not accidentally using Amber python
# (You can comment these 2 lines out once stable)
echo "[pKa filter] python: $(command -v python)"
python -V

# ===============================
# Paths
# ===============================
SCRIPT_DIR="$(cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P)"
RESULTS_DIR="$SCRIPT_DIR/results/lead_generation_output"

# Input/Output
IN="$RESULTS_DIR/generated_smiles_pass_pH6_3.smi"    # input: result filtered at pH 6.4
OUT="$RESULTS_DIR/generated_smiles_pass_pKa.smi"     # final output after pKa filter
LOG="$RESULTS_DIR/pka_log.csv"

mkdir -p "$RESULTS_DIR"

# Avoid "stale output looks like success"
rm -f "$OUT" "$LOG"

# ---- pKa filtering criteria ----
PKA_MIN=4
PKA_MAX=6

# Use absolute path for checkpoint to avoid relative-path surprises
CKPT="$SCRIPT_DIR/datasets/data/pKa/last.ckpt"
# --------------------------------

# Validate inputs early (fail fast, clear error)
if [[ ! -s "$IN" ]]; then
  echo "[pKa filter] ERROR: input file missing/empty: $IN" >&2
  exit 2
fi
if [[ ! -f "$CKPT" ]]; then
  echo "[pKa filter] ERROR: checkpoint not found: $CKPT" >&2
  exit 3
fi

# Run pKa filter
python "$SCRIPT_DIR/utils/predict_pka.py" \
  --in-smi "$IN" \
  --out-csv "$LOG" \
  --out-smi "$OUT" \
  --reg-ckpt "$CKPT" \
  --reg-min "$PKA_MIN" \
  --reg-max "$PKA_MAX"

echo "[pKa filter] Final passed: $(wc -l < "$OUT" | tr -d ' ') molecules -> $OUT"
