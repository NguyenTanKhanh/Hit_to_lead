#!/bin/bash

python vinascreen.py \
  --receptor /home/khanhnt/Desktop/Hit_to_lead/data/Docking/protonation.pdbqt\
  --ligand_dir /home/khanhnt/Desktop/Hit_to_lead/lead_generation_output/pdbqt_files\
  --output_dir docking_results \
  --center 4 50 37 \
  --size 45 45 45 \
  --exhaustiveness 8
