# -*- coding: utf-8 -*-
import os
import argparse
import pandas as pd
from pathlib import Path
from rdkit import Chem
from drugex.training.generators import GraphTransformer
from drugex.data.corpus.vocabulary import VocGraph
from drugex.training.scorers.properties import Property
from drugex.training.scorers.modifiers import ClippedScore
from drugex.training.environment import DrugExEnvironment
from drugex.training.rewards import WeightedSum

def replace_iodine_with_hydrogen(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    rw = Chem.RWMol(mol)
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() == 53:
            atom.SetAtomicNum(1)
    new_mol = rw.GetMol()
    Chem.SanitizeMol(new_mol)
    return Chem.MolToSmiles(new_mol, canonical=True)

def generate_molecules(model_path, frags, num_samples, output_dir, gpu_id=0, min_sa=0.9):
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    use_gpus = (gpu_id,) if gpu_id >= 0 else None

    sa = Property("SA", modifier=ClippedScore(lower_x=7, upper_x=3))
    env = DrugExEnvironment([sa], [1.0], reward_scheme=WeightedSum())

    model = GraphTransformer(voc_trg=VocGraph(), use_gpus=use_gpus)
    model.loadStatesFromFile(model_path)

    print("\n?? Generating molecules...")
    raw = model.generate(input_frags=frags, num_samples=num_samples)

    if isinstance(raw, dict) and 'smiles' in raw:
        smiles = raw['smiles']
    elif isinstance(raw, pd.DataFrame):
        col = next((c for c in raw.columns if c.lower() == 'smiles'), None)
        smiles = raw[col] if col else []
    else:
        smiles = [Chem.MolToSmiles(m) for m in raw]

    smiles = [replace_iodine_with_hydrogen(s) for s in smiles if Chem.MolFromSmiles(s)]

    valid, logs = [], []
    for s in smiles:
        mol = Chem.MolFromSmiles(s)
        if not mol:
            logs.append(f"? Invalid: {s}")
            continue
        score = sa([mol])[0]
        if score > min_sa:
            cano = Chem.MolToSmiles(mol, canonical=True)
            valid.append(cano)
            logs.append(f"? {cano} (SA={score:.3f})")
        else:
            logs.append(f"? {s} (SA={score:.3f})")

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    (Path(output_dir) / "sa_filter_log.txt").write_text("\n".join(logs))
    print("\n".join(logs))

    smi_path = Path(output_dir) / "generated_smiles.smi"
    pd.DataFrame({'SMILES': valid}).to_csv(smi_path, index=False, header=False)
    print(f"\n?? Saved {len(valid)} SMILES to: {smi_path}")
    return str(smi_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--generation_model', required=True)
    parser.add_argument('--frags', nargs='+', required=True)
    parser.add_argument('--num_samples', type=int, default=100)
    parser.add_argument('--output_dir', default='results/lead_generation_output')
    parser.add_argument('--gpu', type=int, default=0)
    parser.add_argument('--min_sa', type=float, default=0.9)
    args = parser.parse_args()

    generate_molecules(
        args.generation_model, args.frags, args.num_samples,
        args.output_dir, args.gpu, args.min_sa
    )
