# -*- coding: utf-8 -*-
import argparse
import os
import sys
import torch
import pandas as pd
from rdkit import Chem

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


def clean_fragments(frag_list, min_sa=0.9):
    """
    Validate and canonicalize SMILES using RDKit.
    Filter fragments with DrugEx SA score > min_sa.
    """
    sa_scorer = Property("SA", modifier=ClippedScore(lower_x=7, upper_x=3))
    valid_frags = []

    for frag in frag_list:
        try:
            mol = Chem.MolFromSmiles(frag)
            if mol:
                sa = sa_scorer([mol])[0]
                if sa > min_sa:
                    canonical = Chem.MolToSmiles(mol, canonical=True)
                    valid_frags.append(canonical)
                    print(f"? Accepted fragment (SA={sa:.3f}): {canonical}")
                else:
                    print(f"?? Skipped fragment (SA={sa:.3f} < {min_sa}): {frag}")
            else:
                print(f"? Invalid SMILES: {frag}")
        except Exception as e:
            print(f"? Error processing '{frag}': {e}")
    return valid_frags


# ---------------------- Argument Parser ----------------------
parser = argparse.ArgumentParser(description='Train DrugEx Graph-based Generator')
parser.add_argument('--gpu', type=int, default=0, help='GPU ID (-1 for CPU)')
parser.add_argument('--frags', nargs='*', help='List of scaffold fragments (SMILES format)')
parser.add_argument('--frag_file', type=str, help='Path to file containing fragment SMILES')
parser.add_argument('--graph_input_folder', type=str, required=True, help='Folder to store encoded data')
parser.add_argument('--batch_size', type=int, default=32, help='Batch size')
parser.add_argument('--epochs', type=int, default=10, help='Number of training epochs')
parser.add_argument('--n_samples', type=int, default=100, help='Samples per scaffold for training')
parser.add_argument('--min_sa', type=float, default=0.9, help='Minimum SA score (0-1) required for fragments')
args = parser.parse_args()


# ---------------------- Combine Fragments ----------------------
all_frags = []

if args.frags:
    all_frags.extend(args.frags)

if args.frag_file:
    if not os.path.exists(args.frag_file):
        print(f"? Fragment file not found: {args.frag_file}")
        sys.exit(1)
    with open(args.frag_file) as f:
        file_frags = [line.strip() for line in f if line.strip()]
        all_frags.extend(file_frags)

if not all_frags:
    print("? No fragments provided. Aborting.")
    sys.exit(1)


# ---------------------- GPU or CPU ----------------------
use_gpus = (args.gpu,) if args.gpu >= 0 else None
os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu) if args.gpu >= 0 else ""

# ---------------------- Reward Environment ----------------------
sascore = Property("SA", modifier=ClippedScore(lower_x=7, upper_x=3))
qed = Property("QED", modifier=ClippedScore(lower_x=0.2, upper_x=0.8))
environment = DrugExEnvironment([sascore, qed], [0.5, 0.5], reward_scheme=WeightedSum())

# ---------------------- Fragment Encoding ----------------------
fragmenter = dummyMolsFromFragments()
encoder = FragmentCorpusEncoder(
    fragmenter=fragmenter,
    encoder=GraphFragmentEncoder(VocGraph(n_frags=4)),
    pairs_splitter=None,
    n_proc=1,
    chunk_size=1
)

# ---------------------- Folders ----------------------
os.makedirs(args.graph_input_folder, exist_ok=True)
training_output_dir = os.path.join("training_graph")
os.makedirs(training_output_dir, exist_ok=True)
encoded_path = os.path.join(args.graph_input_folder, 'scaffolds.tsv')

# ---------------------- Clean Fragments ----------------------
print("?? Cleaning and filtering input fragments...")
cleaned_frags = clean_fragments(all_frags, min_sa=args.min_sa)
if not cleaned_frags:
    print("? No valid fragments after filtering. Aborting.")
    sys.exit(1)

# ---------------------- Apply Encoding ----------------------
dataset = GraphFragDataSet(encoded_path, rewrite=True)
encoder.apply(cleaned_frags, encodingCollectors=[dataset])

# ---------------------- Load Pretrained Models ----------------------
vocab_path = 'data/models/pretrained/Papyrus05.5_graph_trans_PT.vocab'
model_path = 'data/models/pretrained/Papyrus05.5_graph_trans_PT.pkg'

vocabulary = VocGraph.fromFile(vocab_path)

agent = GraphTransformer(voc_trg=vocabulary, use_gpus=use_gpus)
agent.loadStatesFromFile(model_path)

prior = GraphTransformer(voc_trg=vocabulary, use_gpus=use_gpus)
prior.loadStatesFromFile(model_path)

# ---------------------- Explorer ----------------------
explorer = FragGraphExplorer(agent=agent, env=environment, mutate=prior, epsilon=0.1, use_gpus=use_gpus)

# ---------------------- Dataloaders ----------------------
train_loader = GraphFragDataSet(encoded_path).asDataLoader(
    batch_size=args.batch_size, n_samples=args.n_samples)
test_loader = GraphFragDataSet(encoded_path).asDataLoader(
    batch_size=args.batch_size, n_samples=args.n_samples, n_samples_ratio=0.2)

# ---------------------- Training ----------------------
monitor_path = os.path.join(training_output_dir, "scaffolds")
monitor = FileMonitor(monitor_path, save_smiles=True)

explorer.fit(train_loader, test_loader, monitor=monitor, epochs=args.epochs)

# ---------------------- Results ----------------------
model_file = monitor_path + ".pkg"
smiles_file = monitor_path + "_smiles.tsv"

if os.path.isfile(model_file):
    print(f"? Training complete. Model saved to: {model_file}")
else:
    print(f"? Training failed or model not saved.")

if os.path.isfile(smiles_file):
    print(f"?? Sampled SMILES saved to: {smiles_file}")
else:
    print(f"?? No SMILES file found.")
