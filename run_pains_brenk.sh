#!/usr/bin/env bash
set -euo pipefail

# Lấy đường dẫn thư mục của script (dù chạy từ đâu)
SCRIPT_DIR="$(cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P)"
RESULTS_DIR="$SCRIPT_DIR/results/lead_generation_output"

# Input/Output cố định
IN="$RESULTS_DIR/generated_smiles.smi"              # file đầu vào
OUT="$RESULTS_DIR/generated_smiles_pass_pains_brenk.smi"   # file output sau khi lọc PAINS + BRENK

# Tạo thư mục nếu chưa có
mkdir -p "$RESULTS_DIR"

# Kiểm tra file input tồn tại
if [[ ! -f "$IN" ]]; then
  echo "ERROR: Input not found: $IN" >&2
  exit 1
fi

# Gọi Python filter
python utils/pains_and_brenk_filter.py \
  --input "$IN" \
  --output "$OUT"

echo "[PAINS+BRENK] Passed: $(wc -l < "$OUT" | tr -d ' ') molecules -> $OUT"
