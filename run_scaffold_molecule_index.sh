#!/usr/bin/env bash

# Exit on error
set -e

# Define variables
SMILES="NC1=NN(C(=O)C1)C1=CC(C2=NNN=N2)=C(OC2=CC=CC=C2)C=C1"
OUTPUT="results/scaffold/scaffold_indices.png"
WIDTH=500
HEIGHT=500

# Run the Python script
python3 utils/image_scaffold_index.py \
  --smiles "$SMILES" \
  --output "$OUTPUT" \
  --width "$WIDTH" \
  --height "$HEIGHT"

