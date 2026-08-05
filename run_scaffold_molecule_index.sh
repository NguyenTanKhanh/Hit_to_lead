#!/usr/bin/env bash

# Exit on error
set -e

# Define variables
SMILES="CCOC(=O)C4SC3N=c2c(SCCC(=O)O)cc(C(=O)c1ccccc1O)cn2C(=O)C3C4C"
OUTPUT="results/scaffold/scaffold_indices.png"
WIDTH=500
HEIGHT=500

# Run the Python script
python3 utils/image_scaffold_index.py \
  --smiles "$SMILES" \
  --output "$OUTPUT" \
  --width "$WIDTH" \
  --height "$HEIGHT"

