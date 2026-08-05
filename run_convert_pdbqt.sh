#!/bin/bash

# Optional: Activate Conda environment
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate Hit2lead

# === Input file SMILES (bạn có thể đổi ở đây) ===
# Mặc định: ưu tiên active_smiles.smi, nếu không thì dùng generated_smiles.smi
OUTPUT_DIR="results/lead_generation_output"
ACTIVE_SMILES="$OUTPUT_DIR/active_smiles.smi"
ALL_SMILES="$OUTPUT_DIR/generated_smiles.smi"

if [ -f "$ACTIVE_SMILES" ]; then
  SMILES_FILE="$ACTIVE_SMILES"
else
  SMILES_FILE="$ALL_SMILES"
fi

if [ ! -f "$SMILES_FILE" ]; then
  echo "❌ No SMILES file found to convert."
  exit 1
fi

python - <<EOF
from utils.smi_to_pdbqt import convert_smiles_to_pdbqt
convert_smiles_to_pdbqt("$SMILES_FILE", "$OUTPUT_DIR")
EOF

echo "✅ SMILES to PDBQT conversion complete."
