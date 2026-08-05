#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import os
import sys
import csv
import torch
import pandas as pd
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize.rdMolStandardize import LargestFragmentChooser

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


# ---------- Helpers: strip salts/ions and keep the core fragment ----------
remover = SaltRemover()
chooser = LargestFragmentChooser()

def strip_salts_keep_core(smi: str, drop_if_multicomponent: bool = True):
    """
    Parse CXSMILES, remove salts/ions, keep the largest fragment, sanitize it,
    and return the isomeric/canonical SMILES of the core fragment.
    If the processed molecule remains multicomponent and
    drop_if_multicomponent=True, return None.
    If no valid core fragment remains, return None.
    """
    mol = Chem.MolFromSmiles(smi, sanitize=True)
    if not mol:
        return None
    mol = remover.StripMol(mol, dontRemoveEverything=True)
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    mol = chooser.choose(mol)
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    Chem.SanitizeMol(mol)
    core = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    if drop_if_multicomponent and "." in core:
        return None
    return core


def clean_fragments(frag_list, min_sa=0.9):
    """
    Remove salts/ions, keep the largest core fragment, remove any remaining
    multicomponent structures, canonicalize, and filter by SA score.
    """
    sa_scorer = Property("SA", modifier=ClippedScore(lower_x=7, upper_x=3))
    valid_frags = []
    seen = set()

    for frag in frag_list:
        try:
            core = strip_salts_keep_core(frag, drop_if_multicomponent=True)
            if not core:
                print(f"Skipped after salt stripping (empty/ion/multicomponent): {frag}")
                continue

            mol = Chem.MolFromSmiles(core)
            if not mol:
                print(f"Invalid after salt stripping: {frag} -> {core}")
                continue

            sa = sa_scorer([mol])[0]
            if sa <= min_sa:
                print(f"Skipped (SA={sa:.3f} <= {min_sa}): {core}")
                continue

            cano = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            if "." in cano:
                print(f"Skipped (still multicomponent after salt stripping): {cano}")
                continue

            if cano not in seen:
                seen.add(cano)
                valid_frags.append(cano)
                print(f"Accepted (SA={sa:.3f}): {cano}")

        except Exception as e:
            print(f"Error processing '{frag}': {e}")

    return sorted(valid_frags)


# ---------------------- Argument Parser ----------------------
parser = argparse.ArgumentParser(description='Train DrugEx Graph-based Generator')

parser.add_argument('--gpu', type=int, default=0,
                    help='GPU ID (-1 for CPU).')

parser.add_argument('--frags', nargs='*',
                    help='List of scaffold fragments (SMILES format).')

parser.add_argument('--frag_file', type=str,
                    help='Path to file containing fragment SMILES (one per line).')

parser.add_argument('--graph_input_folder', type=str, required=True,
                    help='Folder where `scaffolds.tsv` (encoded data) will be written.')

parser.add_argument('--batch_size', type=int, default=32,
                    help='Batch size for training.')

parser.add_argument('--epochs', type=int, default=10,
                    help='Number of training epochs.')

parser.add_argument('--n_samples', type=int, default=100,
                    help='Samples per scaffold for training.')

parser.add_argument('--min_sa', type=float, default=0.9,
                    help='Minimum SA score (0-1) required to keep a fragment.')

# Allow specifying where to save the final model and SMILES
parser.add_argument('--output_folder', type=str, default='results/training_graph',
                    help='Directory where the final model (.pkg) and sampled SMILES (_smiles.tsv) are saved.')

# Allow specifying the paths to the vocabulary and pretrained model
parser.add_argument('--vocab_path', type=str,
                    default='data/models/pretrained/Papyrus05.5_graph_trans_PT.vocab',
                    help='Path to the pretrained vocabulary file (VocGraph).')

parser.add_argument('--model_path', type=str,
                    default='data/models/pretrained/Papyrus05.5_graph_trans_PT.pkg',
                    help='Path to the pretrained GraphTransformer model (.pkg).')

args = parser.parse_args()


# ---------------------- Combine Fragments ----------------------
all_frags = []

if args.frags:
    all_frags.extend(args.frags)

if args.frag_file:
    if not os.path.exists(args.frag_file):
        print(f"Fragment file not found: {args.frag_file}")
        sys.exit(1)
    with open(args.frag_file) as f:
        file_frags = [line.strip() for line in f if line.strip()]
        all_frags.extend(file_frags)

if not all_frags:
    print("No fragments provided. Aborting.")
    sys.exit(1)


# ---------------------- GPU or CPU ----------------------
if args.gpu >= 0:
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
    use_gpus = (0,)
else:
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
    use_gpus = None


# ---------------------- Reward Environment ----------------------
sascore = Property("SA", modifier=ClippedScore(lower_x=7, upper_x=3))
qed    = Property("QED", modifier=ClippedScore(lower_x=0.2, upper_x=0.8))
environment = DrugExEnvironment(
    [sascore, qed],
    [0.5, 0.5],
    reward_scheme=WeightedSum()
)


# ---------------------- Fragment Encoding ----------------------
fragmenter = dummyMolsFromFragments()
encoder    = FragmentCorpusEncoder(
    fragmenter=fragmenter,
    encoder=GraphFragmentEncoder(VocGraph(n_frags=4)),
    pairs_splitter=None,
    n_proc=1,
    chunk_size=1
)


# ---------------------- Folders ----------------------
# Ensure graph_input_folder exists
os.makedirs(args.graph_input_folder, exist_ok=True)

# Use the user-provided output_folder instead of a hardcoded path
training_output_dir = args.output_folder
os.makedirs(training_output_dir, exist_ok=True)

encoded_path = os.path.join(args.graph_input_folder, 'scaffolds.tsv')


# ---------------------- Clean Fragments ----------------------
print("Cleaning and filtering input fragments...")
cleaned_frags = clean_fragments(all_frags, min_sa=args.min_sa)
if not cleaned_frags:
    print("No valid fragments after filtering. Aborting.")
    sys.exit(1)


# ---------------------- Apply Encoding (safe per-fragment) ----------------------
dataset = GraphFragDataSet(encoded_path, rewrite=True)

skipped = []
ok = 0
for frag in cleaned_frags:
    try:
        encoder.apply([frag], encodingCollectors=[dataset])
        ok += 1
    except Exception as e:
        skipped.append((frag, str(e)))
        print(f"Skipped fragment due to encoding error: {frag}\n   Cause: {e}")

if skipped:
    skip_path = os.path.join(args.graph_input_folder, "skipped_fragments.tsv")
    with open(skip_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["fragment", "error"])
        w.writerows(skipped)
    print(f"==> {len(skipped)} fragment(s) were skipped during encoding. Log: {skip_path}")
else:
    print("==> All fragments encoded successfully.")


# ---------------------- Load Pretrained Models ----------------------
vocab_path  = args.vocab_path
model_path  = args.model_path

if not os.path.exists(vocab_path):
    print(f"Vocabulary file not found: {vocab_path}")
    sys.exit(1)
if not os.path.exists(model_path):
    print(f"Pretrained model file not found: {model_path}")
    sys.exit(1)

vocabulary = VocGraph.fromFile(vocab_path)

agent = GraphTransformer(voc_trg=vocabulary, use_gpus=use_gpus)
agent.loadStatesFromFile(model_path)

prior = GraphTransformer(voc_trg=vocabulary, use_gpus=use_gpus)
prior.loadStatesFromFile(model_path)


# ---------------------- Explorer ----------------------
explorer = FragGraphExplorer(
    agent=agent,
    env=environment,
    mutate=prior,
    epsilon=0.1,
    use_gpus=use_gpus
)


# ---------------------- Dataloaders ----------------------
train_loader = GraphFragDataSet(encoded_path).asDataLoader(
    batch_size=args.batch_size,
    n_samples=args.n_samples
)
test_loader = GraphFragDataSet(encoded_path).asDataLoader(
    batch_size=args.batch_size,
    n_samples=args.n_samples,
    n_samples_ratio=0.2
)


# ---------------------- Training ----------------------
# Monitor files (scaffolds.pkg and scaffolds_smiles.tsv) go under output_folder
monitor_path = os.path.join(training_output_dir, "scaffolds")
monitor = FileMonitor(monitor_path, save_smiles=True)

explorer.fit(train_loader, test_loader, monitor=monitor, epochs=args.epochs)


# ---------------------- Results ----------------------
model_file  = monitor_path + ".pkg"
smiles_file = monitor_path + "_smiles.tsv"

if os.path.isfile(model_file):
    print(f"Training complete. Model saved to: {model_file}")
else:
    print(f"Training failed or model not saved.")

if os.path.isfile(smiles_file):
    print(f"Sampled SMILES saved to: {smiles_file}")
else:
    print(f"No SMILES file found.")