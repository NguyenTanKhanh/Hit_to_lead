

#!/usr/bin/env bash
set -euo pipefail

INPUT_SMI="datasets/FUT8_active_compounds.smi"
PRETRAINED="datasets/data/models/pretrained/Papyrus05.5_graph_trans_PT.pkg"
VOCAB="datasets/data/models/pretrained/Papyrus05.5_graph_trans_PT.vocab"
OUTPUT_DIR="results/training_graph/fine_tune"
EPOCHS=5
BATCH_SIZE=64
GPU=0

mkdir -p "$OUTPUT_DIR"

python utils/fine_tune_drugex_graph.py \
   --input_smi "$INPUT_SMI" \
   --pretrained_model "$PRETRAINED" \
   --vocab_file "$VOCAB" \
   --output_dir "$OUTPUT_DIR" \
   --epochs $EPOCHS \
   --batch_size $BATCH_SIZE \
   --gpu $GPU


### Fine-tune model with a list of compounds intergating some fragments you want to bias the model

##!/usr/bin/env bash
#set -euo pipefail

#INPUT_SMI="datasets/FUT8_active_compounds.smi"
#PRETRAINED="datasets/data/models/pretrained/Papyrus05.5_graph_trans_PT.pkg"
#VOCAB="datasets/data/models/pretrained/Papyrus05.5_graph_trans_PT.vocab"
#OUTPUT_DIR="results/training_graph/fine_tune"
#EPOCHS=5
#BATCH_SIZE=64
#GPU=0
#FRAG_FILE="datasets/my_fragments.smi"

#mkdir -p "$OUTPUT_DIR"

#python utils/fine_tune_drugex_graph.py \
#  --input_smi "$INPUT_SMI" \
#  --pretrained_model "$PRETRAINED" \
#  --vocab_file "$VOCAB" \
#  --output_dir "$OUTPUT_DIR" \
#  --epochs $EPOCHS \
#  --batch_size $BATCH_SIZE \
#  --gpu $GPU \
#  --frag_file "$FRAG_FILE"

