#!/bin/bash

# Optional: Make sure the vina binary is executable
chmod +x utils/vina_1.2.7_linux_x86_64

# Run the Vina docking screen
python utils/vinascreen.py \
  --receptor datasets/data/Docking/protonation.pdbqt \
  --ligand_dir results/lead_generation_output/pdbqt_files \
  --output_dir results/lead_generation_output/docking_results \
  --center 4 50 37 \
  --size 45 45 45 \
  --exhaustiveness 8
