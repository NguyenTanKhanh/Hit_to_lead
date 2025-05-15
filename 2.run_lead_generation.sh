#!/bin/bash

# Optional: activate environment
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate Hit2lead

OUTPUT_DIR="lead_generation_output"

python lead_generation.py \
  --pkg training_graph/scaffolds.pkg \
  --frags \
    "NC1=NN(C(=O)C1)C1=CC(C2=NNN=N2)=C(OC2=CC=CC=C2)C=C1" \
    "CCNC(=O)[N-]S(=O)(=O)C1=CC(=CC=C1OC1=CC=CC=C1)N1N=C(N)CC1=O" \
    "NC1=NN(C(=O)C1)C1=CC=C(OC2=CC=CC=C2)C(=C1)C1=N[N-]N=N1" \
  --num_samples 1000 \
  --output_dir "$OUTPUT_DIR" \
  --gpu 0 \
  --min_sa 0.95

echo "? Lead generation and PDBQT conversion done. Output in: $OUTPUT_DIR/"
