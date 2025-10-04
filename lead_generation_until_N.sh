#!/usr/bin/env bash
set -euo pipefail
set -x

# === Fix environment conflicts ===
# Avoid picking up AMBER's mpi4py and disable Lightning cluster auto-detection
unset PYTHONPATH
unset PYTHONHOME
export PL_DISABLE_ENVIRONMENT_DETECTION=1

# === Config ===
TARGET_N=1000          # s? compound unique c?n thu du?c
BATCH_SIZE=10000         # s? compound generate m?i vòng
SCRIPT_DIR="$(cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P)"
WORK_DIR="$SCRIPT_DIR/results/lead_generation_output"
RESULTS_DIR="$SCRIPT_DIR/results/final_collection"

FINAL_SMI="$RESULTS_DIR/final_pass_${TARGET_N}.smi"
FINAL_CSV="$RESULTS_DIR/final_pass_${TARGET_N}.csv"

mkdir -p "$WORK_DIR" "$RESULTS_DIR"
rm -f "$FINAL_SMI" "$FINAL_CSV"

TOTAL=0
ROUND=0

while [ "$TOTAL" -lt "$TARGET_N" ]; do
  ROUND=$((ROUND+1))
  echo "?? Round $ROUND (Generated so far: $TOTAL / $TARGET_N)"

  # 1. Generate molecules
  python utils/lead_generation.py \
    --generation_model results/training_graph/fine_tune/finetuned.pkg \
    --frags "CC(I)NC(=O)C1=C(I)c2cc(Cl)ccc2OC1I" \
    --num_samples $BATCH_SIZE \
    --output_dir "$WORK_DIR" \
    --gpu 0 \
    --min_sa 0.8

  # 2. PAINS + BRENK filter
  bash "$SCRIPT_DIR/run_pains_brenk.sh" || true
  [[ -s "$WORK_DIR/generated_smiles_pass_pains_brenk.smi" ]] || continue

  # 3. Filters theo pH và ADMET
  bash "$SCRIPT_DIR/run_filter_pH_7.4.sh" || true
  [[ -s "$WORK_DIR/generated_smiles_pass_pH7_4.smi" ]] || continue

  bash "$SCRIPT_DIR/run_filter_pH_6.3.sh" || true
  [[ -s "$WORK_DIR/generated_smiles_pass_pH6_3.smi" ]] || continue

  bash "$SCRIPT_DIR/run_pKa.sh" || true
  [[ -s "$WORK_DIR/generated_smiles_pass_pKa.smi" ]] || continue

  bash "$SCRIPT_DIR/run_predict_activity.sh" || true
  [[ -s "$WORK_DIR/final_pass.smi" ]] || continue

  # 4. Append k?t qu? pass cu?i
  cat "$WORK_DIR/final_pass.smi" >> "$FINAL_SMI"

  # Lo?i b? SMILES trùng l?p
  sort -u "$FINAL_SMI" -o "$FINAL_SMI"

  # CSV cung c?n l?c theo SMILES unique
  if [[ ! -f "$FINAL_CSV" ]]; then
    cp "$WORK_DIR/predict_log.csv" "$FINAL_CSV"
  else
    tail -n +2 "$WORK_DIR/predict_log.csv" >> "$FINAL_CSV"
  fi

  # Gi? CSV d?ng b? v?i SMILES unique
  awk -F, 'NR==FNR {smiles[$1]; next} ($1 in smiles) {print}' "$FINAL_SMI" "$FINAL_CSV" > "$FINAL_CSV.tmp"
  mv "$FINAL_CSV.tmp" "$FINAL_CSV"

  # 5. Update count d?a trên unique SMILES
  TOTAL=$(wc -l < "$FINAL_SMI" | tr -d ' ')
  echo "Current total unique passed: $TOTAL / $TARGET_N"
done

echo "?? Done! Collected $TOTAL unique molecules passing all filters."
echo "?? Final SMILES: $FINAL_SMI"
echo "?? Final CSV:   $FINAL_CSV"
