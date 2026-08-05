#!/usr/bin/env bash
set -euo pipefail
set -x

# === Fix environment conflicts ===
unset PYTHONPATH
unset PYTHONHOME
unset PYTHONUSERBASE
export PYTHONNOUSERSITE=1
export PL_DISABLE_ENVIRONMENT_DETECTION=1

# === Config ===
TARGET_N=10000
BATCH_SIZE=1000
MAX_ROUNDS=1000

SCRIPT_DIR="$(cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P)"
WORK_DIR="$SCRIPT_DIR/results/lead_generation_output"
RESULTS_DIR="$SCRIPT_DIR/results/final_collection"

FINAL_SMI="$RESULTS_DIR/final_pass_${TARGET_N}.smi"

# ? Make sure dirs exist (THIS fixes your error)
mkdir -p "$WORK_DIR" "$RESULTS_DIR"

# Optional: avoid stale output from previous runs
rm -f "$FINAL_SMI"

TOTAL=0
ROUND=0

while [ "$TOTAL" -lt "$TARGET_N" ] && [ "$ROUND" -lt "$MAX_ROUNDS" ]; do
  ROUND=$((ROUND+1))
  echo ">>> Round $ROUND (Generated so far: $TOTAL / $TARGET_N)"

  # 1. Generate molecules
  python utils/lead_generation.py \
    --generation_model results/training_graph/scaffolds.pkg \
    --frags "O=C(C1=C(I)N2C(=O)C3(I)C(I)(N=C2C(SC(I)(I)C(I)(I)C(=O)OI)=C1I)SC(I)(C(=O)OC(I)(I)C(I)(I)I)C3(I)C(I)(I)I)c1ccccc1OI" \
    --num_samples "$BATCH_SIZE" \
    --output_dir "$WORK_DIR" \
    --gpu 0 \
    --min_sa 0.2

  # 2. PAINS + BRENK filter
  bash "$SCRIPT_DIR/run_pains_brenk.sh" || true
  [[ -s "$WORK_DIR/generated_smiles_pass_pains_brenk.smi" ]] || continue

  # 3. Filters theo pH và ADMET
  bash "$SCRIPT_DIR/run_filter_pH_7.4.sh" || true
  [[ -s "$WORK_DIR/generated_smiles_pass_pH7_4.smi" ]] || continue

  bash "$SCRIPT_DIR/run_filter_pH_6.3.sh" || true
  [[ -s "$WORK_DIR/generated_smiles_pass_pH6_3.smi" ]] || continue

  # ? avoid stale pKa output looking like success
  rm -f "$WORK_DIR/generated_smiles_pass_pKa.smi"
  bash "$SCRIPT_DIR/run_pKa.sh" || true
  [[ -s "$WORK_DIR/generated_smiles_pass_pKa.smi" ]] || continue

  # 4. Append result
  # ? defensive: ensure parent dir exists right before append
  mkdir -p "$(dirname "$FINAL_SMI")"
  cat "$WORK_DIR/generated_smiles_pass_pKa.smi" >> "$FINAL_SMI"

  # Remove duplicates
  sort -u "$FINAL_SMI" -o "$FINAL_SMI"

  TOTAL=$(wc -l < "$FINAL_SMI" | tr -d ' ')
  echo "Current total unique passed: $TOTAL / $TARGET_N"
done

if [ "$TOTAL" -ge "$TARGET_N" ]; then
  echo ">>> Done! Reached target: $TOTAL / $TARGET_N"
else
  echo ">>> Stopped: reached MAX_ROUNDS=$MAX_ROUNDS with $TOTAL / $TARGET_N"
fi

echo ">>> Final SMILES: $FINAL_SMI"
