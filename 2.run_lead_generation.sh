#!/bin/bash

# Activate your environment if needed
# source ~/Downloads/yes/envs/DrugEx/bin/activate

# Set default output directory
OUTPUT_DIR="lead_generation_output"

# Run the lead generation and conversion
python lead_generation.py \
  --pkg training_graph/scaffolds.pkg \
  --frags "c1cnccn1" "c1cnccn1.c1ccsc1" \
  --num_samples 100 \
  --output_dir "$OUTPUT_DIR" \
  --gpu 0

echo "Lead generation complete. Results saved to $OUTPUT_DIR/"