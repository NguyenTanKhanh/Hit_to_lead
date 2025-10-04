#!/usr/bin/env bash

# Run DrugEx graph-based training,
# specifying exactly where the pretrained vocab and model live,
# plus where to save our own outputs.

python utils/train_drugex_graph.py \
  --gpu 0 \
  --frag_file "datasets/frags_enamine.smi" \
  --graph_input_folder "datasets/encoded/graph" \
  --batch_size 32 \
  --epochs 5 \
  --n_samples 100 \
  --min_sa 0.8 \
  --vocab_path "datasets/data/models/pretrained/Papyrus05.5_graph_trans_PT.vocab" \
  --model_path "datasets/data/models/pretrained/Papyrus05.5_graph_trans_PT.pkg" \
  --output_folder "results/training_graph"
