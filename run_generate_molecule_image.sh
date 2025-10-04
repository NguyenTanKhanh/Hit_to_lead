#!/usr/bin/env bash

# Exit if any command fails
set -e

# Define paths
SMILES_FILE="results/final_collection/final_pass_1000.smi"
OUTPUT_DIR="results/lead_generation_output/image_lead_molecule"
SCAFFOLD_SMILES="CCNC(=O)C2=Cc1cc(Cl)ccc1OC2"

# Run the Python script
python3 utils/generate_molecule_image.py \
    "$SMILES_FILE" \
    "$OUTPUT_DIR" \
    "$SCAFFOLD_SMILES"
