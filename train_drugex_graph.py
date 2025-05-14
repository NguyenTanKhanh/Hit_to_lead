import argparse
import os
import sys
import torch
import pandas as pd

from drugex.training.scorers.properties import Property
from drugex.training.scorers.modifiers import ClippedScore
from drugex.training.environment import DrugExEnvironment
from drugex.training.rewards import WeightedSum

from drugex.molecules.converters.dummy_molecules import dummyMolsFromFragments
from drugex.data.fragments import FragmentCorpusEncoder, GraphFragmentEncoder
from drugex.data.corpus.vocabulary import VocGraph
from drugex.data.datasets import GraphFragDataSet

from drugex.training.generators import GraphTransformer
from drugex.training.explorers import FragGraphExplorer
from drugex.training.monitors import FileMonitor

# Argument parser
parser = argparse.ArgumentParser(description='Train DrugEx Graph-based Generator')
parser.add_argument('--gpu', type=int, default=0, help='GPU ID (-1 for CPU)')
parser.add_argument('--frags', nargs='+', required=True, help='List of scaffold fragments')
parser.add_argument('--graph_input_folder', type=str, required=True, help='Folder to store encoded data')
parser.add_argument('--batch_size', type=int, default=32, help='Batch size')
parser.add_argument('--epochs', type=int, default=10, help='Number of training epochs')
parser.add_argument('--n_samples', type=int, default=100, help='Samples per scaffold for training')

args = parser.parse_args()

# Set GPU or CPU
use_gpus = (args.gpu,) if args.gpu >= 0 else None
os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu) if args.gpu >= 0 else ""

# Build reward environment
sascore = Property("SA", modifier=ClippedScore(lower_x=7, upper_x=3))
qed = Property("QED", modifier=ClippedScore(lower_x=0.2, upper_x=0.8))
environment = DrugExEnvironment([sascore, qed], [0.5, 0.5], reward_scheme=WeightedSum())

# Fragment encoding
fragmenter = dummyMolsFromFragments()
encoder = FragmentCorpusEncoder(
    fragmenter=fragmenter,
    encoder=GraphFragmentEncoder(VocGraph(n_frags=4)),
    pairs_splitter=None,
    n_proc=1,
    chunk_size=1
)

# Ensure input and output folders exist
os.makedirs(args.graph_input_folder, exist_ok=True)
training_output_dir = os.path.join("training_graph")
os.makedirs(training_output_dir, exist_ok=True)

# Encode fragments
encoded_path = os.path.join(args.graph_input_folder, 'scaffolds.tsv')
dataset = GraphFragDataSet(encoded_path, rewrite=True)
encoder.apply(args.frags, encodingCollectors=[dataset])

# Load vocabulary and pretrained models
vocab_path = '/home/khanhnt/Desktop/Hit_to_lead/data/models/pretrained/graph-trans/Papyrus05.5_graph_trans_PT/Papyrus05.5_graph_trans_PT.vocab'
model_path = '/home/khanhnt/Desktop/Hit_to_lead/data/models/pretrained/graph-trans/Papyrus05.5_graph_trans_PT/Papyrus05.5_graph_trans_PT.pkg'

vocabulary = VocGraph.fromFile(vocab_path)
agent = GraphTransformer(voc_trg=vocabulary, use_gpus=use_gpus)
agent.loadStatesFromFile(model_path)

prior = GraphTransformer(voc_trg=vocabulary, use_gpus=use_gpus)
prior.loadStatesFromFile(model_path)

explorer = FragGraphExplorer(agent=agent, env=environment, mutate=prior, epsilon=0.1, use_gpus=use_gpus)

# DataLoaders
train_loader = GraphFragDataSet(encoded_path).asDataLoader(batch_size=args.batch_size, n_samples=args.n_samples)
test_loader = GraphFragDataSet(encoded_path).asDataLoader(batch_size=args.batch_size, n_samples=args.n_samples, n_samples_ratio=0.2)

# Monitor and Train
monitor_path = os.path.join(training_output_dir, "scaffolds")
monitor = FileMonitor(monitor_path, save_smiles=True)
explorer.fit(train_loader, test_loader, monitor=monitor, epochs=args.epochs)

# Confirm model file exists
model_file = monitor_path + ".pkg"
smiles_file = monitor_path + "_smiles.tsv"
if os.path.isfile(model_file):
    print(f"✅ Training complete. Model saved to: {model_file}")
else:
    print(f"❌ Training failed or model not saved.")

if os.path.isfile(smiles_file):
    print(f"✅ Sampled SMILES saved to: {smiles_file}")
else:
    print(f"⚠️ No SMILES file found.")
