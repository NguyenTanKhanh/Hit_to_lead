#!/bin/bash

# Optional: Make sure the vina binary is executable
chmod +x vina_1.2.7_linux_x86_64

# Run the Vina docking screen
python vinascreen.py \
  --receptor data/Docking/protonation.pdbqt \
  --ligand_dir lead_generation_output/pdbqt_files \
  --output_dir docking_results \
  --center 4 50 37 \
  --size 45 45 45 \
  --exhaustiveness 32
