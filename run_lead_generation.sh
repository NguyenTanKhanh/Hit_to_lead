#!/bin/bash

# Optional: Activate Conda environment
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate Hit2lead

OUTPUT_DIR="results/lead_generation_output"

python utils/lead_generation.py \
  --generation_model results/training_graph/scaffolds.pkg \
  --frags "CC(NC(=O)C1=C(I)c2ccccc2OC1I)C(=O)O" \
  --num_samples 10 \
  --output_dir "$OUTPUT_DIR" \
  --gpu 0 \
  --min_sa 0.8

echo "? Molecule generation complete."
echo "?? Results in: $OUTPUT_DIR/"
