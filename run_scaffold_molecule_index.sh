#!/usr/bin/env bash

# Exit on error
set -e

# Define variables
SMILES="CCNC(=O)C2=Cc1cc(Cl)ccc1OC2"
OUTPUT="results/scaffold/scaffold_indices.png"
WIDTH=500
HEIGHT=500

# Run the Python script
python3 utils/image_scaffold_index.py \
  --smiles "$SMILES" \
  --output "$OUTPUT" \
  --width "$WIDTH" \
  --height "$HEIGHT"

