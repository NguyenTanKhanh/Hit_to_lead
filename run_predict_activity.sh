#!/usr/bin/env bash
set -euo pipefail

# Lấy đường dẫn thư mục của script (dù chạy từ đâu)
SCRIPT_DIR="$(cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P)"
RESULTS_DIR="$SCRIPT_DIR/results/lead_generation_output"

# Input/Output
IN="$RESULTS_DIR/generated_smiles_pass_pKa.smi"    # input: kết quả đã lọc pKa
OUT="$RESULTS_DIR/final_pass.smi"                  # output: kết quả cuối cùng sau activity prediction
LOG="$RESULTS_DIR/predict_log.csv"

# Tạo thư mục nếu chưa có
mkdir -p "$RESULTS_DIR"

# ---- Checkpoint paths & thresholds ----
CLS_CKPT="./datasets/data/classification_model/checkpoint/finetuned_best_model.ckpt"
CLS_THRESHOLD=0.5
REG1_CKPT="./datasets/data/classification_model/checkpoint/first_screening_last.ckpt"
REG1_THRESHOLD=50
REG2_CKPT="./datasets/data/classification_model/checkpoint/second_screening_last.ckpt"
REG2_THRESHOLD=10
# ---------------------------------------

# Chạy activity prediction
python utils/predict_activity.py \
  --in-smi "$IN" \
  --out-csv "$LOG" \
  --out-smi "$OUT" \
  --cls-ckpt "$CLS_CKPT" \
  --cls-threshold "$CLS_THRESHOLD" \
  --reg1-ckpt "$REG1_CKPT" \
  --reg1-threshold "$REG1_THRESHOLD" \
  --reg2-ckpt "$REG2_CKPT" \
  --reg2-threshold "$REG2_THRESHOLD"

echo "[Activity filter] Final passed: $(wc -l < "$OUT" | tr -d ' ') molecules -> $OUT"


# #!/usr/bin/env bash
# set -euo pipefail

# # Lấy đường dẫn thư mục của script (dù chạy từ đâu)
# SCRIPT_DIR="$(cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P)"
# RESULTS_DIR="$SCRIPT_DIR/results/lead_generation_output"

# # Input/Output
# IN="$RESULTS_DIR/generated_smiles_pass_pKa.smi"        # input: kết quả đã lọc pKa
# OUT="$RESULTS_DIR/final_pass_cls_only.smi"             # output: kết quả cuối cùng sau classification
# LOG="$RESULTS_DIR/predict_log_cls_only.csv"

# # Tạo thư mục nếu chưa có
# mkdir -p "$RESULTS_DIR"

# # ---- Checkpoint paths & thresholds ----
# CLS_CKPT="./datasets/data/classification_model/checkpoint/finetuned_best_model.ckpt"
# CLS_THRESHOLD=0.5
# # ---------------------------------------

# # Chạy classification-only
# python utils/predict_activity.py \
#   --in-smi "$IN" \
#   --out-csv "$LOG" \
#   --out-smi "$OUT" \
#   --cls-ckpt "$CLS_CKPT" \
#   --cls-threshold "$CLS_THRESHOLD" \
#   --reg1-ckpt None \
#   --reg2-ckpt None

# echo "[Activity filter - Classification only] Final passed: $(wc -l < "$OUT" | tr -d ' ') molecules -> $OUT"
