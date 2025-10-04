#!/usr/bin/env bash
set -euo pipefail

# Lấy đường dẫn thư mục của script (dù chạy từ đâu)
SCRIPT_DIR="$(cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P)"
RESULTS_DIR="$SCRIPT_DIR/results/lead_generation_output"

# Input/Output cố định
IN="$RESULTS_DIR/generated_smiles_pass_pH7_4.smi"      # input: kết quả đã lọc ở pH 7.4
OUT="$RESULTS_DIR/generated_smiles_pass_pH6_3.smi"     # output: kết quả cuối cùng sau pH 6.3

# Tạo thư mục nếu chưa có
mkdir -p "$RESULTS_DIR"

# ---- Tiêu chí sàng lọc ----
MW_LO=200; MW_HI=800
LOGP_LO=0; LOGP_HI=5
TPSA_LO=0;  TPSA_HI=140
HBD_LO=0;    HBD_HI=5  
HBA_LO=0;    HBA_HI=10
ROTB_LO=0;   ROTB_HI=10
LOGS_LO=-6;  LOGS_HI=0    # logS range (ESOL model)
# ----------------------------

# Chạy lọc
python utils/compound_properties_filter.py \
  --mode filter \
  --in "$IN" \
  --out "$OUT" \
  --ph 6.5 \
  --mw   "$MW_LO" "$MW_HI" \
  --logp "$LOGP_LO" "$LOGP_HI" \
  --tpsa "$TPSA_LO" "$TPSA_HI" \
  --hbd  "$HBD_LO" "$HBD_HI" \
  --hba  "$HBA_LO" "$HBA_HI" \
  --rotb "$ROTB_LO" "$ROTB_HI" \
  --logs "$LOGS_LO" "$LOGS_HI"

echo "[pH 6.3] Final passed: $(wc -l < "$OUT" | tr -d ' ') molecules -> $OUT"
