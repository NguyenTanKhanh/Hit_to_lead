#!/usr/bin/env bash


 set -e

 python3 utils/block_atom.py \
   --smiles "CCOC(=O)C4SC3N=c2c(SCCC(=O)O)cc(C(=O)c1ccccc1O)cn2C(=O)C3C4C" \
   --indices "0,1,5,7,12,13,16,17,27,28,32,33,34" \
   --output-smi results/scaffold/scaffold_block_atom.smi \
   --output-img results/scaffold/scaffold_block_atom.png \
   --width 500 \
   --height 500



