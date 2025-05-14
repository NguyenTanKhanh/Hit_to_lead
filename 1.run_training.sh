#!/bin/bash

# Activate your environment if needed
# source ~/Downloads/yes/envs/DrugEx/bin/activate

# Run the training script with specified arguments
python train_drugex_graph.py \
  --gpu 0 \
  --frags c1cnccn1 c1cnccn1.c1ccsc1 \
  --graph_input_folder datasets/encoded/graph \
  --batch_size 4 \
  --epochs 2 \
  --n_samples 100
