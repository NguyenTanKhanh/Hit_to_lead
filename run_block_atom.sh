#!/usr/bin/env bash


set -e

python3 utils/block_atom.py \
  --smiles "NC1=NN(C(=O)C1)C1=CC(C2=NNN=N2)=C(OC2=CC=CC=C2)C=C1" \
  --indices "21,5,6" \
  --output-smi results/scaffold/scaffold_block_atom.smi \
  --output-img results/scaffold/scaffold_block_atom.png \
  --width 500 \
  --height 500