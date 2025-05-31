#!/bin/bash

# Optional: Activate Conda environment
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate Hit2lead

OUTPUT_DIR="results/lead_generation_output"

python utils/lead_generation.py \
  --generation_model results/training_graph/scaffolds.pkg \
  --prediction_model utils/Classification/checkpoints \
  --frags "NC1=NN(c2c(I)c(I)c(Oc3ccccc3)c(-c3nn[nH]n3)c2I)C(=O)C1" \
  --num_samples 100 \
  --output_dir "$OUTPUT_DIR" \
  --gpu 0 \
  --min_sa 0.8

echo "✅ Lead generation, prediction (if model provided), and PDBQT conversion complete."
echo "📂 Results in: $OUTPUT_DIR/"


