#!/bin/bash

python train_drugex_graph.py \
  --gpu 0 \
  --frag_file "frags_enamine.smi" \
  --graph_input_folder datasets/encoded/graph \
  --batch_size 128 \
  --epochs 100 \
  --n_samples 200 \
  --min_sa 0.8




