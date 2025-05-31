#!/usr/bin/env bash

# Exit if any command fails
set -e

# Define paths
SMILES_FILE="results/lead_generation_output/active_smiles.smi"
OUTPUT_DIR="results/lead_generation_output/image_lead_molecule"
SCAFFOLD_SMILES="NC1=NN(C(=O)C1)C1=CC(C2=NNN=N2)=C(OC2=CC=CC=C2)C=C1"

# Run the Python script
python3 utils/generate_molecule_image.py \
    "$SMILES_FILE" \
    "$OUTPUT_DIR" \
    "$SCAFFOLD_SMILES"
