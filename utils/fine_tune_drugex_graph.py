#!/usr/bin/env python3
import argparse
import os
import shutil
import pandas as pd

from drugex.data.processing import Standardization
from drugex.data.fragments import FragmentCorpusEncoder, GraphFragmentEncoder, FragmentPairsSplitter
from drugex.molecules.converters.fragmenters import Fragmenter
from drugex.data.corpus.vocabulary import VocGraph
from drugex.data.datasets import GraphFragDataSet
from drugex.training.generators import GraphTransformer
from drugex.training.monitors import FileMonitor


def main():
    parser = argparse.ArgumentParser(description="Fine-tune DrugEx GraphTransformer with SMILES input")
    parser.add_argument("--input_smi", required=True, help="Path to .smi file (1 SMILES per line)")
    parser.add_argument("--pretrained_model", required=True, help="Path to pretrained .pkg file")
    parser.add_argument("--vocab_file", required=True, help="Path to pretrained .vocab file")
    parser.add_argument("--output_dir", required=True, help="Directory to save fine-tuned model")
    parser.add_argument("--epochs", type=int, default=3, help="Number of epochs")
    parser.add_argument("--batch_size", type=int, default=64, help="Training batch size")
    parser.add_argument("--gpu", type=int, default=-1, help="GPU id (-1 for CPU)")
    parser.add_argument("--fragments", type=str, nargs="*", default=None,
                        help="Optional list of fragments (SMILES) to constrain generation")
    parser.add_argument("--frag_file", type=str, default=None,
                        help="Optional file containing fragment SMILES (one per line)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # --- GPU / CPU selection ---
    use_gpus = (args.gpu,) if args.gpu >= 0 else None
    if args.gpu >= 0:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
        print(f"✅ Using GPU {args.gpu}")
    else:
        print("✅ Using CPU")

    # --- Load SMILES ---
    with open(args.input_smi) as f:
        smiles = [line.strip().split()[0] for line in f if line.strip()]
    print(f"Loaded {len(smiles)} SMILES from {args.input_smi}")

    # --- Standardize ---
    standardizer = Standardization(n_proc=4)
    smiles = standardizer.apply(smiles)

    # --- Prepare encoder ---
    encoder = FragmentCorpusEncoder(
        fragmenter=Fragmenter(4, 4, 'brics'),
        encoder=GraphFragmentEncoder(VocGraph(n_frags=4)),
        pairs_splitter=FragmentPairsSplitter(0.1, 100),  # 10% test split
        n_proc=4
    )

    train_path = os.path.join(args.output_dir, "train.tsv")
    test_path = os.path.join(args.output_dir, "test.tsv")

    train = GraphFragDataSet(train_path, rewrite=True)
    test = GraphFragDataSet(test_path, rewrite=True)

    encoder.apply(list(smiles), encodingCollectors=[test, train])

    # --- Load pretrained model ---
    vocab = VocGraph.fromFile(args.vocab_file)
    pretrained = GraphTransformer(voc_trg=vocab, use_gpus=use_gpus)
    pretrained.loadStatesFromFile(args.pretrained_model)

    # --- Monitor ---
    monitor = FileMonitor(
        os.path.join(args.output_dir, "finetuned"),
        save_smiles=True,
        reset_directory=True
    )

    # --- Train ---
    pretrained.fit(
        train.asDataLoader(args.batch_size),
        test.asDataLoader(args.batch_size),
        epochs=args.epochs,
        monitor=monitor
    )

    # --- Save vocab (copy) ---
    vocab_out = os.path.join(args.output_dir, "finetuned.vocab")
    shutil.copy(args.vocab_file, vocab_out)
    print(f"✅ Saved vocab to {vocab_out}")

    # --- Optional: generate with fragment constraints ---
    frag_list = []
    if args.fragments:
        frag_list.extend(args.fragments)
    if args.frag_file and os.path.exists(args.frag_file):
        with open(args.frag_file) as f:
            frag_list.extend([line.strip() for line in f if line.strip()])

    if frag_list:
        print(f"✅ Generating with {len(frag_list)} input fragments...")
        generated = pretrained.generate(input_frags=frag_list, num_samples=50, batch_size=args.batch_size)
        frag_out = os.path.join(args.output_dir, "finetuned_fragments.tsv")
        pd.DataFrame(generated).to_csv(frag_out, sep="\t", index=False)
        print(f"✅ Saved fragment-constrained generation to {frag_out}")

    print("🎉 Fine-tuning complete.")


if __name__ == "__main__":
    main()
