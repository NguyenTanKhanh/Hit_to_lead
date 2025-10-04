#!/usr/bin/env bash


 set -e

 python3 utils/block_atom.py \
   --smiles "CCNC(=O)C2=Cc1cc(Cl)ccc1OC2" \
   --indices "1,6,15" \
   --output-smi results/scaffold/scaffold_block_atom.smi \
   --output-img results/scaffold/scaffold_block_atom.png \
   --width 500 \
   --height 500




